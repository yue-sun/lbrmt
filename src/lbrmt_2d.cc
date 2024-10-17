// #include <algorithm>
#include <cstring>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
// #include <vector>

#include "flow_model.hh"
#include "lbrmt_2d.hh"
#include "omp.h"
#include "convert_unit.hh"

/** Construct to the lbrmt class.
 * \param[in] rho_f_ the fluid density.
 * \param[in] Re_ the Reynolds number.
 * \param[in] tau_ the relaxation time.
 * \param[in] epsilon_ the blur zone half-width.
 * \param[in] nx_ the number of grid points in the x direction.
 * \param[in] ny_ the number of grid points in the y direction.
 * \param[in] ax_ the lower bound of simulation domain in the x direction.
 * \param[in] bx_ the upper bound of simulation domain in the x direction.
 * \param[in] ay_ the lower bound of simulation domain in the y direction.
 * \param[in] by_ the upper bound of simulation domain in the y direction.
 * \param[in] bctype_ the boundary condition type.
 * \param[in] nobjs_ the number of all solid objects.
 * \param[in] outdir_ the output directory.
 * \param[in] outdir_ the output directory for multimaps.
 * \param[in] obsfile_ an optional binary file of boolean flags denoting rigid obstacles.
 * \param[in] ux0_ the initial fluid velocity in the x direction.
 * \param[in] uy0_ the initial fluid velocity in the y direction. */
lbrmt_2d::lbrmt_2d(double rho_f_,double Re_,double tau_,double epsilon_,
    int nx_,int ny_,double ax_,double bx_,double ay_,double by_,
    int bctype_,int nobjs_,char *outdir_, char *outdir_mmap_,
    const char *obsfile_,double ux0_,double uy0_) :
    rho_f(rho_f_), Re(Re_), tau(tau_), nu((tau-0.5)/3.), epsilon(epsilon_),
    nx(nx_), nxb(nx+4), ny(ny_), nyb(ny+4), ax(ax_), bx(bx_), ay(ay_), by(by_),
    bctype(bctype_), ux0(ux0_), uy0(uy0_), nobjs(nobjs_), time(0.), counter_ov(0), 
    fbase(new lattice[nxb*nyb]), f(fbase+2*nxb+2), obs(new bool[nxb*nyb]), of_list(new obj_field*[nobjs]) {
    // Set up output directory
    set_dir(outdir_,outdir_mmap_);
    // Set up rigid obstacles boolean array
    if(obsfile_!=NULL) {
        FILE *fb=safe_fopen(obsfile_,"rb");
        safe_fread(obs,sizeof(bool),nxb*nyb,fb,"obstacles");
        fclose(fb);
    }
}

/** Destory the class objects. */
lbrmt_2d::~lbrmt_2d() {
    // Free object-related data
    for(int i=nobjs-1;i>=0;i--) delete of_list[i];
    delete [] of_list;
    delete [] obs;
    // Free ref_map objects
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;
            if(fbase[k].max_mmap>0) delete [] fbase[k].mmap;
        }
    }
    // Free the main grid data structure
    delete [] fbase;
    // Free output buffers
    delete [] outdir;
    delete [] outdir_mmap;
}

/** Add an object to the simulation, creating a new obj_field class for all solid fields.
 * \param[in] o the pointer to the object to add.
 * \param[in] mc the material constants for the solid object.
 * \param[in] oid the object ID. */
void lbrmt_2d::add_object(object *o,mat_const *mc,int oid) {
    // Create the new obj_field class on the list
    of_list[oid]=new obj_field(o,mc,epsilon);
}

/** Initialize reference map for solid nodes. */
void lbrmt_2d::init_remap() {
#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;

            // Initialize mmap to default values
            // Make resolve_overflow() call consistent with empty mmap list
            fbase[k].counter_mmap=0;
            fbase[k].max_mmap=0;

            // Loop oer all solid objects, q: object ID
            for(int q=0;q<nobjs;q++) {
                object *o=of_list[q]->obj;
                // Initialize the obj_field and objects
                o->start(*this,of_list[q]);
                // Initialize the reference map fields to be the identity map
                double x=ax+i-1.5;
                double y=ay+j-1.5;                
                // Initialize the level set values
                double XX,YY;
                o->pre_stress_setup(0.);
                o->transform(x,y,XX,YY); // (x,y) world space, (X,Y) object space
                double phi=o->phi(XX,YY);
                // Create ref_map if within the blur zone
                if(phi<(epsilon<2.?(epsilon<1.?8.*epsilon:4.*epsilon):2.*epsilon)) {
                    // Initially allocate length of 2 for the mmap array
                    // Since assumed all solids are not overlapping, only 1 ref_map is present
                    if(fbase[k].counter_mmap<1) fbase[k].mmap=new ref_map[4];
                    // Store world space reference map
                    fbase[k].mmap[fbase[k].counter_mmap]=ref_map(x,y,phi,q,epsilon);
                    // Increment counter of mmap and set the maximum array length
                    fbase[k].counter_mmap+=1;
                    fbase[k].max_mmap=4;
                }
            }
        }
    }
}

/** Initialize hydro variables of all nodes.
 * \param[in] rho the macroscopic density.
 * \param[in] macro_ux the x component of macroscopic velocity.
 * \param[in] macro_uy the y component of macroscopic velocity.
 * \param[in] fx_fluid the x component of external force density on the fluid nodes.
 * \param[in] fy_fluid the y component of external force density on the fluid nodes. */
void lbrmt_2d::init_hydro(double rho,double macro_ux,double macro_uy,double fx_fluid,double fy_fluid) {
#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;

            // Initialize macroscopic quantities to be of fluid
            double rho_h=rho_f,ux_h=macro_ux,uy_h=macro_uy,fx_h=fx_fluid,fy_h=fy_fluid;
            double drho=0.,sum_lambda=0.;
            // Add transition to solids by looping over each object in mmap on one node
            for(int q=0;q<fbase[k].counter_mmap;q++) {
                ref_map &rmap=fbase[k].mmap[q];  
                // Initialize the solid fraction and get solid density from the corresponding object ID's obj_field
                // HACK: this only works in initialization b/c solids are not overlapping
                rmap.solid_fraction();
                double rho_s=of_list[rmap.id]->rho_s;
                // Compute the Heaviside blur density
                drho+=rmap.lambda*(rho_s-rho_f);
                // Compute the Heaviside force density
                double fx_solid=0.,fy_solid=0.;
                double x=ax+i-1.5;
                double y=ay+j-1.5;
                // Transform xx and yy to the correct object space
                object *o=of_list[rmap.id]->obj;
                o->accel(x,y,rmap.X,rmap.Y,rmap.phi,fx_solid,fy_solid);
                fx_h=rmap.lambda*fx_solid+fx_fluid;
                fy_h=rmap.lambda*fy_solid+fy_fluid;
            }
            if(sum_lambda>1.) drho/=sum_lambda;
            rho_h+=drho;

            // Open channel B.C.
            if(bctype==2) {
                // Inlet: set velocity
                if(i<2) ux_h=Re*nu/ny;
                // Outlet: set density
                if(i>nx-1) rho_h=rho_f;
            }

            // Convergence study: Taylor–Green vortex decay
            if(bctype==8) {
                double x=ax+i-1.5;
                double y=ay+j-1.5;
                int scale=ny/40;
                taylor_green(0.,x,y,scale,rho_h,ux_h,uy_h);
            }
            // Initialize hydro variables of each node
            fbase[k].init_hydro(rho_h,ux_h,uy_h,fx_h,fy_h);

            // Add obstacle flags for rigid obstacles
            if(obs[k]) {
                fbase[k].obs=true;
                fbase[k].ux=0.;
                fbase[k].uy=0.;
            } else fbase[k].obs=false;
        }
    }    
}

/** Initialize populations of all nodes. */
void lbrmt_2d::init_pop() {
#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;
            // Call local routine at each node
            fbase[k].init_pop();
        }
    }
}

/** Calculate updated hydro variables. */
void lbrmt_2d::hydro() {
#pragma omp parallel for
    for(int j=2;j<nyb-2;j++) {
        for(int i=2;i<nxb-2;i++) {
            int k=j*nxb+i;
            // Call local routine at each node
            fbase[k].hydro();
        }
    }
}

/** Calculate updated reference map field. */
void lbrmt_2d::remap() {
    // Clean up ref_map objects
#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;
            for(int q=0;q<fbase[k].counter_mmap;) {
                int max_ll=ceil(1.5+sqrt(2.)*epsilon+7);
                if(fbase[k].mmap[q].cc>max_ll-1) {
                    fbase[k].mmap[q]=fbase[k].mmap[--fbase[k].counter_mmap];
                } else q++;
            }
        }
    }

#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;

            // Loop over each object in mmap for one node
            for(int q=0;q<fbase[k].counter_mmap;q++) {
                ref_map &rmap=fbase[k].mmap[q];
                int id=rmap.id;
                object *o=of_list[id]->obj;
                double xx,yy;
                o->transform(rmap.X,rmap.Y,xx,yy); // (rmap.X,rmap.Y) world space, (xx,yy) object space
                double phi=o->phi(xx,yy);
                // Update the level set value
                rmap.phi=phi;
                rmap.cc=rmap.phi<0.?0:-1;
                // Update solid fraction lambda of each object
                rmap.solid_fraction();
                // Reset extrapolation flag for the next extrapolation
                rmap.extrap_reset=0;
            }

            // Reset the components of solid stress at each node
            // It should be in compute_stress(), moved here b/c it needs to cover the entire domain
            fbase[k].sla=0.;fbase[k].slb=0.;fbase[k].slc=0.;fbase[k].sld=0.;
            fbase[k].sba=0.;fbase[k].sbb=0.;fbase[k].sbc=0.;fbase[k].sbd=0.;
            fbase[k].Fx=0.;fbase[k].Fy=0.;fbase[k].total_lambda=0.;
            fbase[k].invdetF_left=0.;fbase[k].invdetF_bottom=0.;fbase[k].avg_invdetF=100000001.;
        }
    }

#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;
            double u=f[k].ux;
            double v=f[k].uy;
            double Xx,Xy,Yx,Yy;

            // Loop over each object in mmap for one node
            for(int q=0;q<f[k].counter_mmap;q++) {
                if(f[k].mmap[q].phi<epsilon) {
                    int id=f[k].mmap[q].id;
                    // Second-order upwinding ENO method
                    if(u>0.) {
                        int id_right=match_id(i+1,j,k+1,id),id_left=match_id(i-1,j,k-1,id),id_left_2=match_id(i-2,j,k-2,id);
                        if(id_right==nobjs||id_left==nobjs||id_left_2==nobjs) {
                            fprintf(stderr,"Can't match in remap() u>0, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
                            exit(1);
                        }
                        Xx=eno2(f[k+1].mmap[id_right].X,f[k].mmap[q].X,f[k-1].mmap[id_left].X,f[k-2].mmap[id_left_2].X);
                        Yx=eno2(f[k+1].mmap[id_right].Y,f[k].mmap[q].Y,f[k-1].mmap[id_left].Y,f[k-2].mmap[id_left_2].Y);
                    } else {
                        int id_left=match_id(i-1,j,k-1,id),id_right=match_id(i+1,j,k+1,id),id_right_2=match_id(i+2,j,k+2,id);
                        if(id_left==nobjs||id_right==nobjs||id_right_2==nobjs) {
                            fprintf(stderr,"Can't match in remap() u<0, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
                            exit(1);
                        }
                        Xx=-eno2(f[k-1].mmap[id_left].X,f[k].mmap[q].X,f[k+1].mmap[id_right].X,f[k+2].mmap[id_right_2].X);
                        Yx=-eno2(f[k-1].mmap[id_left].Y,f[k].mmap[q].Y,f[k+1].mmap[id_right].Y,f[k+2].mmap[id_right_2].Y);
                    }
                    if(v>0.) {
                        int id_top=match_id(i,j+1,k+nxb,id),id_bottom=match_id(i,j-1,k-nxb,id),id_bottom_2=match_id(i,j-2,k-2*nxb,id);
                        if(id_top==nobjs||id_bottom==nobjs||id_bottom_2==nobjs) {
                            fprintf(stderr,"Can't match in remap() v>0, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
                            exit(1);
                        }
                        Xy=eno2(f[k+nxb].mmap[id_top].X,f[k].mmap[q].X,f[k-nxb].mmap[id_bottom].X,f[k-2*nxb].mmap[id_bottom_2].X);
                        Yy=eno2(f[k+nxb].mmap[id_top].Y,f[k].mmap[q].Y,f[k-nxb].mmap[id_bottom].Y,f[k-2*nxb].mmap[id_bottom_2].Y);
                    } else {
                        int id_bottom=match_id(i,j-1,k-nxb,id),id_top=match_id(i,j+1,k+nxb,id),id_top_2=match_id(i,j+2,k+2*nxb,id);
                        if(id_bottom==nobjs||id_top==nobjs||id_top_2==nobjs) {
                            fprintf(stderr,"Can't match in remap() v<0, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
                            exit(1);
                        }
                        Xy=-eno2(f[k-nxb].mmap[id_bottom].X,f[k].mmap[q].X,f[k+nxb].mmap[id_top].X,f[k+2*nxb].mmap[id_top_2].X);
                        Yy=-eno2(f[k-nxb].mmap[id_bottom].Y,f[k].mmap[q].Y,f[k+nxb].mmap[id_top].Y,f[k+2*nxb].mmap[id_top_2].Y);
                    }
                    if(f[k].mmap[q].phi<0.) {
                        f[k].mmap[q].cX=u*Xx+v*Xy;  
                        f[k].mmap[q].cY=u*Yx+v*Yy;
                    } else {
                        f[k].mmap[q].cX=0.;  
                        f[k].mmap[q].cY=0.;
                    }
                    
                    // Compute determinant of deformation gradient tensor
                    double detF=calc_detF(Xx,Xy,Yx,Yy);
                    f[k].mmap[q].detF=detF;
                } else {
                    f[k].mmap[q].cX=0.;  
                    f[k].mmap[q].cY=0.;
                    f[k].mmap[q].detF=1.;
                }
            }
        }
    }

#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;
            for(int q=0;q<fbase[k].counter_mmap;q++) {
                ref_map &rmap=fbase[k].mmap[q];
                // Call local routine at each node
                rmap.update();
            }
        }
    }
}

/** Extrapolate the reference map values in the blur zone positive direction. */
void lbrmt_2d::extrap() {
    // Loop over the entire grid to find the first layer and store to each object
#pragma omp parallel for
    for(int j=2;j<nyb-2;j++) {
        for(int i=2;i<nxb-2;i++) {
            int k=j*nxb+i;
            // Loop over existing mmap array of a node
            for(int q=0;q<fbase[k].counter_mmap;q++) {
                ref_map &rmap=fbase[k].mmap[q];
                int id=rmap.id;
                if(rmap.phi<0.) {
                    std::vector<int> &l=of_list[id]->l;
                    build_layer(i-1,j,k-1,1,l,id);
                    build_layer(i,j-1,k-nxb,1,l,id);
                    build_layer(i+1,j,k+1,1,l,id);
                    build_layer(i,j+1,k+nxb,1,l,id);
                }
            }
        }
    }

    // Loop over all objects to set extrapolated reference map values for the first layer
#pragma omp parallel for
    for(int q=0;q<nobjs;q++) {
        // Retrieve the list of extrapolation of each object
        std::vector<int> &l=of_list[q]->l;
        int p_max=l.size()/2;
        for(int p=0;p<p_max;p++) {
            int k=l[2*p],id=l[2*p+1],i=k%nxb,j=(k-i)/nxb;
            set_layer(i,j,k,id);
        }
    }

    // Loop over all objects
    // TODO: loop over the entire grid if objects are few
    for(int q=0;q<nobjs;q++) {
        // Retrieve the list of extrapolation of each object
        std::vector<int> &l=of_list[q]->l,&l2=of_list[q]->l2;

        // Build subsequent layers based on existing layers
        int max_ll=ceil(1.5+sqrt(2.)*epsilon+7);
        for(int ll=2;ll<=max_ll;ll++) {
            std::vector<int> &al=ll&1?l2:l,&bl=ll&1?l:l2;
            int p_max=al.size()/2;
       
            // Loop over the previous layer to build the current one
#pragma omp parallel for
            for(int p=0;p<p_max;p++) {
                int k=al[2*p],id=al[2*p+1],i=k%nxb,j=(k-i)/nxb;
                if(i>0) build_layer(i-1,j,k-1,ll,bl,id);
                if(j>0) build_layer(i,j-1,k-nxb,ll,bl,id);
                if(i<nxb-1) build_layer(i+1,j,k+1,ll,bl,id);
                if(j<nyb-1) build_layer(i,j+1,k+nxb,ll,bl,id);
            }

            // Set extrapolated reference map values for the current layer
            p_max=bl.size()/2;
#pragma omp parallel for
            for(int p=0;p<p_max;p++) {
                int k=bl[2*p],id=bl[2*p+1],i=k%nxb,j=(k-i)/nxb;
                set_layer(i,j,k,id);
            }
            al.clear();
        }

        // Loop over all objects to reset the extrapolation layer vectors
        of_list[q]->l.clear();
        of_list[q]->l2.clear();
    }
}

/** Calculate the total solid fractions for all objects. */
void lbrmt_2d::total_solid_fraction() {
#pragma omp parallel for
    for(int j=-2;j<ny+2;j++) {
        for(int i=-2;i<nx+2;i++) {
            int k=j*nxb+i;

            double temp_tsf=0.;
            // Loop over each object in mmap for one node
            for(int q=0;q<f[k].counter_mmap;q++) {
                ref_map &rmap=f[k].mmap[q];
                if(rmap.phi<epsilon) temp_tsf+=rmap.lambda;
            }
            f[k].total_lambda=temp_tsf;
        }
    }
}

/** Calculate all stress for all solid nodes. */
void lbrmt_2d::compute_stress() {
    // Perform any object-specific calculations prior to the stress computation
    for(int q=0;q<nobjs;q++) {
        object *o=of_list[q]->obj;
        o->pre_stress_setup(time);
    }

#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;

            // Loop over each object in mmap for one node
            for(int q=0;q<f[k].counter_mmap;q++) {
                ref_map &rmap_center=f[k].mmap[q];

                // Construct edge stresses 1 grid point larger than blur zone
                if(rmap_center.phi<epsilon) {
                    int id=rmap_center.id;
                    object *o=of_list[id]->obj;
                    double G=of_list[id]->G;

                    // Set left edge solid stress
                    solid_stress<true>(i,j,k,id,rmap_center,o,G);
                    // Set down edge solid stress
                    solid_stress<false>(i,j,k,id,rmap_center,o,G);

                    // Collision detection is triggered only when mmap has more than 2 ref_map
                    if(f[k].counter_mmap>1) {
                        // Compute left edge collision stress
                        collision_stress<true>(i,j,k,id,q,rmap_center,o,G);
                        // Compute down edge collision stress
                        collision_stress<false>(i,j,k,id,q,rmap_center,o,G);
                    }
                }
            }
        }
    }
}

/** Helper function to calculate penalized-incompressible solid stress tensor.
 * \param[in] i the x index.
 * \param[in] j the y index.
 * \param[in] k the index of solid node in the f array.
 * \param[in] id the object ID.
 * \param[in] rmap_center the reference to the center ref_map object.
 * \param[in] o the pointer to the object.
 * \param[in] G the shear modulus of the object. */
template<bool left>
void lbrmt_2d::solid_stress(int i,int j,int k,int id,ref_map &rmap_center,object *o,double G) {
    double phiv;                // the half-edge level set value
    double X,Y,Xx,Xy,Yx,Yy;     // the components of the Jacobian
    double BP11,BP12,BP21,BP22; // the components of the deviatoric part of FFT
    double invdetF;
    double J;                   // the determinant of the deformation gradient F and pow(J,-5./3)
    double sf;                  // the coefficient to smooth out blur zone solid stress at half-edge

    // Compute the edge-based reference map and deformation gradient tensor 
    if(left) {
        // Instantiate the five neighboring nodes for stress computation (six in total)
        ref_map &rmap_left=f[k-1].mmap[match_id(i-1,j,k-1,id)];
        ref_map &rmap_lower_left=f[k-1-nxb].mmap[match_id(i-1,j-1,k-1-nxb,id)];
        ref_map &rmap_upper_left=f[k-1+nxb].mmap[match_id(i-1,j+1,k-1+nxb,id)];
        ref_map &rmap_lower_center=f[k-nxb].mmap[match_id(i,j-1,k-nxb,id)];
        ref_map &rmap_upper_center=f[k+nxb].mmap[match_id(i,j+1,k+nxb,id)];
        if(match_id(i-1,j,k-1,id)==nobjs||match_id(i-1,j-1,k-1-nxb,id)==nobjs||match_id(i-1,j+1,k-1+nxb,id)==nobjs
            ||match_id(i,j-1,k-nxb,id)==nobjs||match_id(i,j+1,k+nxb,id)==nobjs) {
            fprintf(stderr,"Can't match in solid_stress() left, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
            exit(1);
        }
        // Calculate half-edge reference map
        phiv=0.5*(rmap_center.phi+rmap_left.phi);
        X=0.5*(rmap_center.X+rmap_left.X);
        Y=0.5*(rmap_center.Y+rmap_left.Y);
        // Calculate the components of the Jacobian using second-order finite difference method
        Xx=rmap_center.X-rmap_left.X;
        Yx=rmap_center.Y-rmap_left.Y;
        Xy=(rmap_upper_center.X+rmap_upper_left.X-rmap_lower_center.X-rmap_lower_left.X)*0.25;
        Yy=(rmap_upper_center.Y+rmap_upper_left.Y-rmap_lower_center.Y-rmap_lower_left.Y)*0.25;
    } else {
        // Instantiate the five neighboring nodes for stress computation (six in total)
        ref_map &rmap_lower_center=f[k-nxb].mmap[match_id(i,j-1,k-nxb,id)];
        ref_map &rmap_left=f[k-1].mmap[match_id(i-1,j,k-1,id)];
        ref_map &rmap_right=f[k+1].mmap[match_id(i+1,j,k+1,id)];
        ref_map &rmap_lower_left=f[k-1-nxb].mmap[match_id(i-1,j-1,k-1-nxb,id)];
        ref_map &rmap_lower_right=f[k+1-nxb].mmap[match_id(i+1,j-1,k+1-nxb,id)];
        if(match_id(i,j-1,k-nxb,id)==nobjs||match_id(i-1,j,k-1,id)==nobjs||match_id(i+1,j,k+1,id)==nobjs
            ||match_id(i-1,j-1,k-1-nxb,id)==nobjs||match_id(i+1,j-1,k+1-nxb,id)==nobjs) {
            fprintf(stderr,"Can't match in solid_stress() right, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
            exit(1);
        }
        // Calculate half-edge reference map
        phiv=0.5*(rmap_center.phi+rmap_lower_center.phi);
        X=0.5*(rmap_center.X+rmap_lower_center.X);
        Y=0.5*(rmap_center.Y+rmap_lower_center.Y);
        // Calculate the components of the Jacobian using second-order finite difference method
        Xx=(rmap_right.X+rmap_lower_right.X-rmap_left.X-rmap_lower_left.X)*0.25;
        Yx=(rmap_right.Y+rmap_lower_right.Y-rmap_left.Y-rmap_lower_left.Y)*0.25;
        Xy=rmap_center.X-rmap_lower_center.X;
        Yy=rmap_center.Y-rmap_lower_center.Y;
    }
        
    // Calculate the determinant of F and B'
    // Also apply an active component F_a from the object's actuate() function
    calc_fft(X,Y,Xx,Xy,Yx,Yy,J,BP11,BP12,BP21,BP22,o);
    // Compute the local the smoothing coefficient
    sf=heaviside(phiv);
    // Compute the inverse of determinant of F
    invdetF=calc_invdetF(Xx,Xy,Yx,Yy);

    if(left) {
        // Compute the smoothed solid stress left components
        f[k].sla+=sf*(G*BP11);
        f[k].slb+=sf*(G*BP12);
        f[k].slc+=sf*(G*BP21);
        f[k].sld+=sf*(G*BP22);
        f[k].invdetF_left+=sf*invdetF;
    } else {
        // Compute the smoothed solid stress bottom components
        f[k].sba+=sf*(G*BP11);
        f[k].sbb+=sf*(G*BP12);
        f[k].sbc+=sf*(G*BP21);
        f[k].sbd+=sf*(G*BP22);
        f[k].invdetF_bottom+=sf*invdetF;
    }
}

/** Helper function to calculate collision stress based on contact forces between solid objects.
 * \param[in] i the x index.
 * \param[in] j the y index.
 * \param[in] k the index of solid node in the f array.
 * \param[in] id the object ID.
 * \param[in] q the index of the mmap array.
 * \param[in] rmap_center the reference to the center ref_map object.
 * \param[in] o the pointer to the object.
 * \param[in] G the shear modulus of the object. */
template<bool left>
void lbrmt_2d::collision_stress(int i,int j,int k,int id,int q,ref_map &rmap_center,object *o,double G) {
    double p1,p1x,p1y;
    // Compute the edge-based reference map and deformation gradient tensor 
    if(left) {
        // Instantiate the five neighboring nodes for stress computation (six in total)
        ref_map &rmap_left=f[k-1].mmap[match_id(i-1,j,k-1,id)];
        ref_map &rmap_lower_left=f[k-1-nxb].mmap[match_id(i-1,j-1,k-1-nxb,id)];
        ref_map &rmap_upper_left=f[k-1+nxb].mmap[match_id(i-1,j+1,k-1+nxb,id)];
        ref_map &rmap_lower_center=f[k-nxb].mmap[match_id(i,j-1,k-nxb,id)];
        ref_map &rmap_upper_center=f[k+nxb].mmap[match_id(i,j+1,k+nxb,id)];
        if(match_id(i-1,j,k-1,id)==nobjs||match_id(i-1,j-1,k-1-nxb,id)==nobjs||match_id(i-1,j+1,k-1+nxb,id)==nobjs
            ||match_id(i,j-1,k-nxb,id)==nobjs||match_id(i,j+1,k+nxb,id)==nobjs) {
            fprintf(stderr,"Can't match in collision_stress() left, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
            exit(1);
        }
        // Compute the edge-based level set and its gradient using second-order finite difference method
        p1=0.5*(rmap_center.phi+rmap_left.phi);
        p1x=rmap_center.phi-rmap_left.phi;
        p1y=(rmap_upper_center.phi+rmap_upper_left.phi-rmap_lower_center.phi-rmap_lower_left.phi)*0.25;
    } else {
        // Instantiate the five neighboring nodes for stress computation (six in total)
        ref_map &rmap_lower_center=f[k-nxb].mmap[match_id(i,j-1,k-nxb,id)];
        ref_map &rmap_left=f[k-1].mmap[match_id(i-1,j,k-1,id)];
        ref_map &rmap_right=f[k+1].mmap[match_id(i+1,j,k+1,id)];
        ref_map &rmap_lower_left=f[k-1-nxb].mmap[match_id(i-1,j-1,k-1-nxb,id)];
        ref_map &rmap_lower_right=f[k+1-nxb].mmap[match_id(i+1,j-1,k+1-nxb,id)];
        if(match_id(i,j-1,k-nxb,id)==nobjs||match_id(i-1,j,k-1,id)==nobjs||match_id(i+1,j,k+1,id)==nobjs
            ||match_id(i-1,j-1,k-1-nxb,id)==nobjs||match_id(i+1,j-1,k+1-nxb,id)==nobjs) {
            fprintf(stderr,"Can't match in collision_stress() bottom, t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
            exit(1);
        }
        // Compute the edge-based level set and its gradient using second-order finite difference method
        p1=0.5*(rmap_center.phi+rmap_lower_center.phi);
        p1x=(rmap_right.phi+rmap_lower_right.phi-rmap_left.phi-rmap_lower_left.phi)*0.25;
        p1y=rmap_center.phi-rmap_lower_center.phi;
    }

    // Loop over all other objects, considering pairs (i,j) where i<j
    for(int q2=q+1;q2<f[k].counter_mmap;q2++) {
        double p2,p2x,p2y;
        int id2=f[k].mmap[q2].id;
        ref_map rmap_center2=f[k].mmap[q2];
        if((rmap_center2.phi<epsilon)) {
            if(left) {
                // Instantiate the five neighboring nodes for stress computation (six in total)
                ref_map &rmap_left2=f[k-1].mmap[match_id(i-1,j,k-1,id2)];
                ref_map &rmap_lower_left2=f[k-1-nxb].mmap[match_id(i-1,j-1,k-1-nxb,id2)];
                ref_map &rmap_upper_left2=f[k-1+nxb].mmap[match_id(i-1,j+1,k-1+nxb,id2)];
                ref_map &rmap_lower_center2=f[k-nxb].mmap[match_id(i,j-1,k-nxb,id2)];
                ref_map &rmap_upper_center2=f[k+nxb].mmap[match_id(i,j+1,k+nxb,id2)];
                if(match_id(i-1,j,k-1,id2)==nobjs||match_id(i-1,j-1,k-1-nxb,id2)==nobjs||match_id(i-1,j+1,k-1+nxb,id2)==nobjs
                    ||match_id(i,j-1,k-nxb,id2)==nobjs||match_id(i,j+1,k+nxb,id2)==nobjs) {
                    fprintf(stderr,"Can't match in collision_stress() other left, t=%g, id=%d, i=%d, j=%d\n",time,id2,i,j);
                    exit(1);
                }
                // Compute the edge-based level set and its gradient using second-order finite difference method
                p2=0.5*(rmap_center2.phi+rmap_left2.phi);
                p2x=rmap_center2.phi-rmap_left2.phi;
                p2y=(rmap_upper_center2.phi+rmap_upper_left2.phi-rmap_lower_center2.phi-rmap_lower_left2.phi)*0.25;
            } else {
                // Instantiate the five neighboring nodes for stress computation (six in total)
                ref_map &rmap_lower_center2=f[k-nxb].mmap[match_id(i,j-1,k-nxb,id2)];
                ref_map &rmap_left2=f[k-1].mmap[match_id(i-1,j,k-1,id2)];
                ref_map &rmap_right2=f[k+1].mmap[match_id(i+1,j,k+1,id2)];
                ref_map &rmap_lower_left2=f[k-1-nxb].mmap[match_id(i-1,j-1,k-1-nxb,id2)];
                ref_map &rmap_lower_right2=f[k+1-nxb].mmap[match_id(i+1,j-1,k+1-nxb,id2)];
                if(match_id(i,j-1,k-nxb,id2)==nobjs||match_id(i-1,j,k-1,id2)==nobjs||match_id(i+1,j,k+1,id2)==nobjs
                    ||match_id(i-1,j-1,k-1-nxb,id2)==nobjs||match_id(i+1,j-1,k+1-nxb,id2)==nobjs) {
                    fprintf(stderr,"Can't match in collision_stress() other bottom, t=%g, id=%d, i=%d, j=%d\n",time,id2,i,j);
                    exit(1);
                }
                // Compute the edge-based level set and its gradient using second-order finite difference method
                p2=0.5*(rmap_center2.phi+rmap_lower_center2.phi);
                p2x=(rmap_right2.phi+rmap_lower_right2.phi-rmap_left2.phi-rmap_lower_left2.phi)*0.25;
                p2y=rmap_center2.phi-rmap_lower_center2.phi;
            }
        
            // Compute collision stress
            double norx,nory,nn,fac;
            double G2=of_list[id2]->G; // G of the other object, G1 is passed as an input
            double eta=8.0;            // Dimensionless constant

            norx=p2x-p1x;
            nory=p2y-p1y;
            nn=norx*norx+nory*nory;
            fac=eta*fmin(contact_force(p1),contact_force(p2))*(G+G2);
            fac=nn<1e-12?0.:fac/nn;
            
            if(left) {
                // Stress left
                f[k].sla-=fac*0.5*(norx*norx-nory*nory);
                f[k].slb-=fac*norx*nory;
                f[k].slc-=fac*norx*nory;
                f[k].sld-=fac*0.5*(nory*nory-norx*norx);
            } else {
                // Stress bottom
                f[k].sba-=fac*0.5*(norx*norx-nory*nory);
                f[k].sbb-=fac*norx*nory;
                f[k].sbc-=fac*norx*nory;
                f[k].sbd-=fac*0.5*(nory*nory-norx*norx);
            }
        }
    }
}

/** Calculate the equilibrium populations. */
void lbrmt_2d::equilibrium() {
#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;

            // Compute target density difference
            double drho=0.,sum_rho_s=0.,sum_detF=0.;
            // Loop over each object in mmap for one node
            for(int q=0;q<fbase[k].counter_mmap;q++) {
                ref_map &rmap=fbase[k].mmap[q];
                if(rmap.phi<epsilon) {
                    double rho_s=of_list[rmap.id]->rho_s;
                    sum_rho_s+=rmap.lambda*rho_s;    
                    sum_detF+=rmap.detF;                
                }
            }
            double lambda=fbase[k].total_lambda;
            if(lambda>1.) {
                sum_rho_s/=lambda;
                sum_detF/=lambda;
            }
            if(lambda>0.) {
                double rho_t=sum_rho_s+(1.-lambda)*rho_f;
                drho=rho_t-rho_f;
            }
            
            // Call local routine at each node
            fbase[k].equilibrium(drho);
        }
    }
}

/** Collide and calculate the populations for the next timestep. */
void lbrmt_2d::collide() {
    double invtau=1./tau;
#pragma omp parallel for
    for(int j=-2;j<ny+2;j++) {
        for(int i=-2;i<nx+2;i++) {
            int k=j*nxb+i;
            // Call local routine at each node
            f[k].collide(invtau);
        }
    }
}

/** Add the external force density to account for forces acting on the fluid nodes (e.g. gravity, pressure gradient)
 * and the force density for forces acting on the solid nodes (e.g. external force, solid stress).
 * \param[in] fx_fluid the x component of external force density on the fluid nodes.
 * \param[in] fy_fluid the y component of external force density on the fluid nodes.
 * \param[in] bc_type the type of boundary conditions. */
void lbrmt_2d::force_density(double fx_fluid,double fy_fluid,int bctype) {
#pragma omp parallel for
    for(int j=-2;j<ny+2;j++) {
        for(int i=-2;i<nx+2;i++) {
            int k=j*nxb+i;

            double fx_solid=0.,fy_solid=0.;
            double fx_fluid_mask=0.,fy_fluid_mask=0.;
            f[k].ss=0.;

            // Loop over each object in mmap for one node
            for(int q=0;q<f[k].counter_mmap;q++) {
                ref_map &rmap=f[k].mmap[q];
                if(rmap.phi<epsilon) {
                    double x=ax+i+0.5;
                    double y=ay+j+0.5;                    
                    // Apply extra accelerations to the object (e.g. anchor force, gravity)
                    double accx=0.,accy=0.;
                    object *o=of_list[rmap.id]->obj;
                    o->accel(x,y,rmap.X,rmap.Y,rmap.phi,accx,accy);
                    fx_solid+=rmap.lambda*accx;
                    fy_solid+=rmap.lambda*accy;
                    // Mask out fluid force
                    fx_fluid_mask+=rmap.lambda*fx_fluid;
                    fy_fluid_mask+=rmap.lambda*fy_fluid;

                    if(bctype==3) {
                        double force_scale=50.0*of_list[rmap.id]->G; // wall repulsive force
                        // Strength of the wall force function, w:support range;
                        const double w=epsilon<=2.?5.:2.*epsilon;
                        const double winv=1/w,k_rep=force_scale*winv;
                        double phiv=rmap.phi/epsilon;
                        // Left and right walls
                        if(x<=ax+w) {
                            fx_solid+=k_rep*delta_func_wall((x-ax)*winv)
                                    *(phiv>-1.?(phiv>1.?0:0.5*(1-phiv*1.)):1);
                        } else if(x>=bx-w) {
                            fx_solid-=k_rep*delta_func_wall((bx-x)*winv)
                                    *(phiv>-1.?(phiv>1.?0:0.5*(1-phiv*1.)):1);
                        }
                        // Bottom and top walls
                        if(y<=ay+w) {
                            fy_solid+=k_rep*delta_func_wall((y-ay)*winv)
                                    *(phiv>-1.?(phiv>1.?0:0.5*(1-phiv*1.)):1);
                        } else if(y>=by-w) {
                            fy_solid-=k_rep*delta_func_wall((by-y)*winv)
                                    *(phiv>-1.?(phiv>1.?0:0.5*(1-phiv*1.)):1);
                        }
                    }

                    // Store magnitude of solid stress (von Mises)
                    if(rmap.phi<epsilon&&(i<nx+1&&j<ny+1)) {
                        double s11=0.25*(f[k].sla+f[k+1].sla+f[k].sba+f[k+nxb].sba);
                        double s12=0.25*(f[k].slb+f[k+1].slb+f[k].sbb+f[k+nxb].sbb);
                        double s21=0.25*(f[k].slc+f[k+1].slc+f[k].sbc+f[k+nxb].sbc);
                        double s22=0.25*(f[k].sld+f[k+1].sld+f[k].sbd+f[k+nxb].sbd);
                        double svm=sqrt(1.5*(s11*s11+s22*s22+s12*s12+s21*s21));
                        f[k].ss+=rmap.lambda*svm;
                    } 
                }
            }

            // If the solid fraction is greater than 1, then scale the force down by this fraction
            if(f[k].total_lambda>1.) {
                fx_solid/=f[k].total_lambda;
                fy_solid/=f[k].total_lambda;
                fx_fluid_mask/=f[k].total_lambda;
                fy_fluid_mask/=f[k].total_lambda;
            }

            // Add divergence of solid stress and inverse of detF
            double divsx=0.,divsy=0.,invdetF=0.;
            if(i<nx+1 && j<ny+1) {
                divsx=f[k+1].sla-f[k].sla+f[k+nxb].sbb-f[k].sbb;
                divsy=f[k+1].slc-f[k].slc+f[k+nxb].sbd-f[k].sbd;
                invdetF=0.25*(f[k+1].invdetF_left+f[k].invdetF_left
                             +f[k+nxb].invdetF_bottom+f[k].invdetF_bottom);
            }
            if(f[k].total_lambda>1.) {
                divsx/=f[k].total_lambda;
                divsy/=f[k].total_lambda;
                invdetF/=f[k].total_lambda;
            }

            // Set net macroscopic force density
            f[k].Fx=fx_fluid-fx_fluid_mask+fx_solid+divsx;
            f[k].Fy=fy_fluid-fy_fluid_mask+fy_solid+divsy;
            f[k].avg_invdetF=invdetF;
        }
    }
}

/** Add forcing term. */
void lbrmt_2d::force() {
#pragma omp parallel for
    for(int j=0;j<nyb;j++) {
        for(int i=0;i<nxb;i++) {
            int k=j*nxb+i;
            // Call local routine at each node
            fbase[k].force(tau);
        }
    }
}

/** Set the boundary conditions.
 * \param[in] bc_type the type of boundary conditions.*/
void lbrmt_2d::bc(int bctype) {
    switch(bctype) {
        // Periodic channel (force-driven)
        case 0: pbc_lbm();break;
        // Open channel (pressure)
        case 2: obc_lbm();break;
        // Fully no-slip box
        case 3: fnbc_lbm();break;
        // Fully periodic box
        case 4: fpbc_lbm();break;
        // Lid-driven cavity
        case 7: ldc_lbm();break;
        // Taylor–Green vortex
        case 8: fpbc_lbm();break;
    }
}

/** Stream post-collision populations along D2Q9 directions.
 *  Temporary storage of streamed populations in f_eq. */
void lbrmt_2d::stream() {
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;
            // Standard streaming
            f[k].feq0=f[k].f0;
            f[k].feq1=f[k-1].f1;
            f[k].feq2=f[k-nxb].f2;
            f[k].feq3=f[k+1].f3;
            f[k].feq4=f[k+nxb].f4;
            f[k].feq5=f[k-nxb-1].f5;
            f[k].feq6=f[k-nxb+1].f6;
            f[k].feq7=f[k+nxb+1].f7;
            f[k].feq8=f[k+nxb-1].f8;
        }
    }
#pragma omp parallel for
    for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) f[j*nxb+i].update_f();
}

/** Add obstacles. */
void lbrmt_2d::obstacle() {
    // Impose no-slip one grid next to the obstacle boundary
    // No need to behave like bc(), i.e. this function happens after streaming
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;

            if(f[k+1].obs) f[k].feq3=f[k].f1;
            if(f[k+nxb].obs) f[k].feq4=f[k].f2;
            if(f[k-1].obs) f[k].feq1=f[k].f3;
            if(f[k-nxb].obs) f[k].feq2=f[k].f4;
            if(f[k+1+nxb].obs) f[k].feq7=f[k].f5;
            if(f[k-1+nxb].obs) f[k].feq8=f[k].f6;
            if(f[k-1-nxb].obs) f[k].feq5=f[k].f7;
            if(f[k+1-nxb].obs) f[k].feq6=f[k].f8;

            if(f[k].obs) {
                f[k].ux=0.;f[k].uy=0.;
            }
        }
    }
#pragma omp parallel for
    for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) f[j*nxb+i].update_f();
}

/** Initialization function called in the main routine.
 * \param[in] rho_f fluid density.
 * \param[in] macro_ux the x component of macroscopic velocity.
 * \param[in] macro_uy the y component of macroscopic velocity.
 * \param[in] bc_type the type of boundary conditions.
 * \param[in] fx_fluid the x component of external force density on the fluid nodes.
 * \param[in] fy_fluid the y component of external force density on the fluid nodes. */
void lbrmt_2d::initialize(double rho_f,double macro_ux,double macro_uy,int bctype,double fx_fluid,double fy_fluid) {
    init_remap();
    init_hydro(rho_f,macro_ux,macro_uy,fx_fluid,fy_fluid);
    equilibrium();
    init_pop();
    bc(bctype);
    obstacle();
}

/** Main solver to step forward in time using the 2D LBRMT.
 * \param[in] niters the number of iterations.
 * \param[in] nout the number of equally-spaced output frames.
 * \param[in] bc_type the type of boundary conditions. 
 * \param[in] fx_fluid the x component of external force density on the fluid nodes.
 * \param[in] fy_fluid the y component of external force density on the fluid nodes. */
void lbrmt_2d::solve(int niters,int nout,int bctype,double fx_fluid,double fy_fluid) {int skip=niters/(nout-1);
    int fr=0;
    output(fr++);
    for(int t=1;t<niters+1;t++) {
        // Compute external force density
        remap();
        extrap();
        total_solid_fraction();
        compute_stress();
        force_density(fx_fluid,fy_fluid,bctype);
        // Update hydro moments
        equilibrium();
        force();
        collide();
        bc(bctype);
        stream();
        obstacle();
        hydro();
        // diag_total_rho();
        if(t%skip==0) output(fr++);
        time+=1.0;
    }
}

/** Main solver to step forward in time using the 2D LBRMT with timing.
 * \param[in] niters the number of iterations.
 * \param[in] nout the number of equally-spaced output frames.
 * \param[in] bc_type the type of boundary conditions. 
 * \param[in] fx_fluid the x component of external force density on the fluid nodes.
 * \param[in] fy_fluid the y component of external force density on the fluid nodes. */
void lbrmt_2d::solve_time(int niters,int nout,int bctype,double fx_fluid,double fy_fluid) {
    // Define timing variables
    double time_bc=0.,time_stream=0.,time_hydro=0.,time_remap=0.,time_extrap=0.,
           time_tsf=0.,time_stress=0.,time_eq=0.,time_collide=0.,time_force=0.;
    clock_t cpu_begin;

    int skip=niters/(nout-1);
    int fr=0;
    output(fr++);
    for(int t=1;t<niters+1;t++) {
        cpu_begin=clock();
        bc(bctype);
        time_bc+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        cpu_begin=clock();
        stream();
        time_stream+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        time+=1.0;

        cpu_begin=clock();
        hydro();
        time_hydro+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        if(t%skip==0) output(fr++);

        cpu_begin=clock();
        remap();
        time_remap+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        cpu_begin=clock();
        extrap();
        time_extrap+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        // diag_count_ref_map();

        cpu_begin=clock();
        total_solid_fraction();
        time_tsf+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        cpu_begin=clock();
        compute_stress();
        time_stress+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        cpu_begin=clock();
        equilibrium();
        time_eq+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        cpu_begin=clock();
        collide();
        time_collide+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);

        cpu_begin=clock();
        force_density(fx_fluid,fy_fluid,bctype);
        time_force+=double(clock()-cpu_begin)/double(CLOCKS_PER_SEC);
    }
    double time_rmt=time_remap+time_extrap+time_tsf+time_stress;
    double time_lbm=time_bc+time_stream+time_hydro+time_eq+time_collide+time_force;
    double time_total=time_rmt+time_lbm;
    printf("==========================================================================\n"
           "C++ routine                             Percentage     CPU time\n"
           "==========================================================================\n"
           "bc()                                    %.8f     %g\n"
           "stream()                                %.8f     %g\n"
           "hydro()                                 %.8f     %g\n"
           "remap()                                 %.8f     %g\n"
           "extrap()                                %.8f     %g\n"
           "total_solid_fraction()                  %.8f     %g\n"
           "compute_stress()                        %.8f     %g\n"
           "equilibrium()                           %.8f     %g\n"
           "collide()                               %.8f     %g\n"
           "force()                                 %.8f     %g\n"
           "==========================================================================\n"
           "RMT routines                            %.8f     %g\n"
           "LBM routines                            %.8f     %g\n"
           "ALL routines                                           %g\n"
           "==========================================================================\n",
           time_bc/time_total,time_bc,time_stream/time_total,time_stream,time_hydro/time_total,time_hydro,
           time_remap/time_total,time_remap,time_extrap/time_total,time_extrap,time_tsf/time_total,time_tsf,
           time_stress/time_total,time_stress,time_eq/time_total,time_eq,time_collide/time_total,time_collide,
           time_force/time_total,time_force,time_rmt/time_total,time_rmt,time_lbm/time_total,time_lbm,time_total);
}

/** Periodic channel boundary conditions for lbm variables.
 *  West periodic, east periodic.
 *  North no-slip, south no-slip.  */
void lbrmt_2d::pbc_lbm() {
    // West inlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb-1;
        f[k].f1=f[k+nx].f1;
        f[k].f5=f[k+nx].f5;
        f[k].f8=f[k+nx].f8;
    }
    // East outlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb+nx;
        f[k].f3=f[k-nx].f3;
        f[k].f6=f[k-nx].f6;
        f[k].f7=f[k-nx].f7;
    }
    // North solid
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=nxb*ny+i;
        f[k].f4=f[k-nxb].f2;
        f[k].f7=f[k-nxb-1].f5;
        f[k].f8=f[k-nxb+1].f6;
    }
    // South solid
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=i-nxb;
        f[k].f2=f[k+nxb].f4;
        f[k].f5=f[k+nxb+1].f7;
        f[k].f6=f[k+nxb-1].f8;
    }
    // Northwest corner bounce-back
    f[nxb*ny-1].f8=f[nxb*(ny-1)].f6;
    // Northeast corner bounce-back
    f[nxb*ny+nx].f7=f[nxb*(ny-1)+nx-1].f5;
    // Southwest corner bounce-back
    f[-1*nxb-1].f5=f[0].f7;
    // Southeast corner bounce-back
    f[-1*nxb+nx].f6=f[nx-1].f8;
}

/** Open channel boundary conditions for lbm variables.
 *  West inlet, east outlet.
 *  North no-slip, south no-slip.  */
void lbrmt_2d::obc_lbm() {
    // West inlet
#pragma omp parallel for
    for(int j=-1;j<ny+1;j++) {
        int k=j*nxb-1;
        double ux=Re*nu/ny,uy=0.;
        double uxsq=ux*ux,uysq=uy*uy,usq=uxsq+uysq;
        double rho=f[k+1].rho;
        f[k].f1=w1*rho*(1.+ux/cs2+0.5*uxsq/cs4-0.5*usq/cs2);
        f[k].f5=w2*rho*(1.+( ux+uy)/cs2+0.5*( ux+uy)*( ux+uy)/cs4-0.5*usq/cs2);
        f[k].f8=w2*rho*(1.+( ux-uy)/cs2+0.5*( ux-uy)*( ux-uy)/cs4-0.5*usq/cs2);
    }
    // East outlet
#pragma omp parallel for
    for(int j=-1;j<ny+1;j++) {
        int k=j*nxb+nx;
        double ux=f[k-1].ux,uy=f[k-1].uy;
        double uxsq=ux*ux,uysq=uy*uy,usq=uxsq+uysq;
        double rho=rho_f;
        f[k].f3=w1*rho*(1.-ux/cs2+0.5*uxsq/cs4-0.5*usq/cs2);
        f[k].f6=w2*rho*(1.+(-ux+uy)/cs2+0.5*(-ux+uy)*(-ux+uy)/cs4-0.5*usq/cs2);
        f[k].f7=w2*rho*(1.+(-ux-uy)/cs2+0.5*(-ux-uy)*(-ux-uy)/cs4-0.5*usq/cs2);
    }
    // North solid
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=nxb*ny+i;
        f[k].f4=f[k-nxb].f2;
        f[k].f7=f[k-nxb-1].f5;
        f[k].f8=f[k-nxb+1].f6;
    }
    // South solid
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=i-nxb;
        f[k].f2=f[k+nxb].f4;
        f[k].f5=f[k+nxb+1].f7;
        f[k].f6=f[k+nxb-1].f8;
    }
}

/** Fully no-slip box boundary conditions for lbm variables.
 *  West no-slip, east no-slip.
 *  North no-slip, south no-slip.  */
void lbrmt_2d::fnbc_lbm() {
    // West inlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb-1;
        f[k].f1=f[k+1].f3;
        f[k].f5=f[k+nxb+1].f7;
        f[k].f8=f[k-nxb+1].f6;
    }
    // East outlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb+nx;
        f[k].f3=f[k-1].f1;
        f[k].f6=f[k+nxb-1].f8;
        f[k].f7=f[k-nxb-1].f5;
    }
    // North solid
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=nxb*ny+i;
        f[k].f4=f[k-nxb].f2;
        f[k].f7=f[k-nxb-1].f5;
        f[k].f8=f[k-nxb+1].f6;
    }
    // South solid
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=i-nxb;
        f[k].f2=f[k+nxb].f4;
        f[k].f5=f[k+nxb+1].f7;
        f[k].f6=f[k+nxb-1].f8;
    }
    // Northwest corner bounce-back
    f[nxb*ny-1].f8=f[nxb*(ny-1)].f6;
    // Northeast corner bounce-back
    f[nxb*ny+nx].f7=f[nxb*(ny-1)+nx-1].f5;
    // Southwest corner bounce-back
    f[-1*nxb-1].f5=f[0].f7;
    // Southeast corner bounce-back
    f[-1*nxb+nx].f6=f[nx-1].f8;
}

/** Fully periodic boundary condition.
 *  West periodic, east periodic.
 *  North periodic, south periodic.  */
void lbrmt_2d::fpbc_lbm(){
    // West inlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb-1;
        f[k].f1=f[k+nx].f1;
        f[k].f5=f[k+nx].f5;
        f[k].f8=f[k+nx].f8;
    }
    // East outlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb+nx;
        f[k].f3=f[k-nx].f3;
        f[k].f6=f[k-nx].f6;
        f[k].f7=f[k-nx].f7;
    }
    // North outlet
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=nxb*ny+i;
        int kk=i;
        f[k].f4=f[kk].f4;
        f[k].f7=f[kk].f7;
        f[k].f8=f[kk].f8;
    }
    // South outlet
#pragma omp parallel for
    for(int i=0;i<nx;i++) {
        int k=i-nxb;
        int kk=nxb*(ny-1)+i;
        f[k].f2=f[kk].f2;
        f[k].f5=f[kk].f5;
        f[k].f6=f[kk].f6;
    }
    // Northwest corner periodic
    f[nxb*ny-1].f8=f[nx-1].f8;
    // Northeast corner periodic
    f[nxb*ny+nx].f7=f[0].f7;
    // Southwest corner periodic
    f[-1*nxb-1].f5=f[nxb*(ny-1)+nx-1].f5;
    // Southeast corner periodic
    f[-1*nxb+nx].f6=f[nxb*(ny-1)].f6;
}

/** Lid-driven cavity boundary conditions for lbm variables.
 *  West no-slip, east no-slip.
 *  North prescribed velocity, south no-slip. */
void lbrmt_2d::ldc_lbm() {
    // Wall boundary velocity
    double ux_wall=Re/ny*nu;
    // West inlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb-1;
        f[k].f1=f[k+1].f3;
        f[k].f5=f[k+nxb+1].f7;
        f[k].f8=f[k-nxb+1].f6;
    }
    // East outlet
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        int k=j*nxb+nx;
        f[k].f3=f[k-1].f1;
        f[k].f6=f[k+nxb-1].f8;
        f[k].f7=f[k-nxb-1].f5;
    }
    // North solid
#pragma omp parallel for
    for(int i=-1;i<nx+1;i++) {
        int k=nxb*ny+i;
        // Moving wall bounce-back, p.p 180 Kruger's book
        f[k].f4=f[k-nxb].f2;
        f[k].f7=f[k-nxb-1].f5-1./6*rho_f*ux_wall;
        f[k].f8=f[k-nxb+1].f6+1./6*rho_f*ux_wall;
    }
    // South solid
#pragma omp parallel for
    for(int i=-1;i<nx+1;i++) {
        int k=i-nxb;
        f[k].f2=f[k+nxb].f4;
        f[k].f5=f[k+nxb+1].f7;
        f[k].f6=f[k+nxb-1].f8;
    }
}

/** Set output directory.
 * \param[in] outdir_ the output directory name.
 * \param[in] outdir_mmap_ the output directory name for multimaps. */
void lbrmt_2d::set_dir(char *outdir_,char *outdir_mmap_) {
    // Output for standard physical fields
    size_t l=strlen(outdir_)+1;
    outdir=new char[2*l+32];
    memcpy(outdir,outdir_,sizeof(char)*l);
    outbuf=outdir+l;
    // Create output directory, if it doesn't exist.
    mkdir(outdir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Output for multimaps
    size_t l_mmap=strlen(outdir_mmap_)+1;
    outdir_mmap=new char[2*l_mmap+32];
    memcpy(outdir_mmap,outdir_mmap_,sizeof(char)*l_mmap);
    outbuf_mmap=outdir_mmap+l_mmap;
    // Create output directory, if it doesn't exist.
    mkdir(outdir_mmap,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
}

/** Write to file.
 * \param[in] fr the frame number */
void lbrmt_2d::output(int fr) {
    // Header line prints the timestep and system dimensions
    sprintf(outbuf,"%s/fr_%04d.txt",outdir,fr);
    fp=safe_fopen(outbuf,"w");
    fprintf(fp,"0 0 0 0 0 %d %d %d\n",fr,nx,ny);
    // Header line prints the timestep and system dimensions for multimaps
    sprintf(outbuf_mmap,"%s/fr_%04d.txt",outdir_mmap,fr);
    fp_mmap=safe_fopen(outbuf_mmap,"w");
    fprintf(fp_mmap,"0 0 0 0 %d %d %d %d\n",fr,nx,ny,nobjs);
    // Output the physical fields
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;
            fprintf(fp,"%d %d %g %g %g %g %g %g\n",i,j,f[k].rho,f[k].ux,f[k].uy,f[k].Fx,f[k].Fy,f[k].ss);
            for(int q=0;q<f[k].counter_mmap;q++) {
                ref_map &rmap=f[k].mmap[q];
                fprintf(fp_mmap,"%d %d %d %g %g %g %d %g\n",rmap.id,i,j,rmap.X,rmap.Y,rmap.phi,rmap.cc,rmap.detF); 
            }
        }
    }
    fclose(fp);
    fclose(fp_mmap);
    printf("Output %d done.\n",fr);
}

/** Second-order ENO upwinding method. */
inline double lbrmt_2d::eno2(double v0,double v1,double v2,double v3) {
    return fabs(v0-2*v1+v2)<fabs(v1-2*v2+v3)?0.5*(v0-v2):1.5*v1-2*v2+0.5*v3;
}

inline void lbrmt_2d::build_layer(int i,int j,int k,int ll,std::vector<int> &l,int id) {
#pragma omp critical
    {
    int mkid=match_base_id(k,id);
    bool mmap_exist=mkid!=nobjs;

    // Case 1: mmap exists and matching
    if(mmap_exist) {
        bool is_fluid=fbase[k].mmap[mkid].cc==-1;
        bool extrap_rest=fbase[k].mmap[mkid].extrap_reset==0;
        // Prevent double extrapolating
        if(is_fluid&&extrap_rest) {
            fbase[k].mmap[mkid].cc=ll;
            // Add to current list, turn on reset flag
            l.push_back(k);
            l.push_back(id);
            fbase[k].mmap[mkid].extrap_reset++;
        }
    // Case 2: no mmap exists or matching
    } else {
        int d=fbase[k].counter_mmap++;
        // Case 2.1: max_length can store one more ref_map
        if(d+1<fbase[k].max_mmap) {
            // Create a new ref_map at the end
            fbase[k].mmap[d]=ref_map(1000000.,1000000.,10000.,id,epsilon);
        // Case 2.2: max_length cannot store one more ref_map
        } else {
            // Allocate more memory by doubling max_mmap
            int omax_mmap=fbase[k].max_mmap;            // Old max_mmap
            int nmax_mmap=(omax_mmap>0)?omax_mmap<<1:4; // New (doubled) max_mmap
            // Create a new array twice the length of the old one
            ref_map* nmmap=new ref_map[nmax_mmap];
            // Copy the old array to the new one
            memcpy((void*)nmmap,(void*)fbase[k].mmap,sizeof(ref_map)*omax_mmap);
            nmmap[d]=ref_map(1000000.,1000000.,10000.,id,epsilon);
            delete [] fbase[k].mmap;
            fbase[k].mmap=nmmap;
            fbase[k].max_mmap=nmax_mmap;
        }
        fbase[k].mmap[d].cc=ll;
        // Add to current list, turn on reset flag
        l.push_back(k);
        l.push_back(id);
        fbase[k].mmap[d].extrap_reset++;
    }
    }
}

/** Resolve the overflow buffer of the extrapolation zone.
 * \param[in] ll the extrapolation layer ID.
 * \param[in] l the vector that stores all (k,id) pairs for nodes in the current layer. */
inline void lbrmt_2d::resolve_overflow(int ll,std::vector<int> &l) {
    int p_max=ov_buffer.size()/2;
    printf("overflow length: %d\n",p_max);
    for(int p=0;p<p_max;p++) {
        int k=ov_buffer[2*p],id=ov_buffer[2*p+1];
        printf("k: %d, id: %d, ll: %d\n",k,id,ll);

        // Prevent double-counting for nodes have been added to the overflow buffer multiple times
        int mkid=match_base_id(k,id);
        if(mkid==nobjs)
        {
            // Allocate more memory by doubling max_mmap
            int omax_mmap=fbase[k].max_mmap;            // Old max_mmap
            int nmax_mmap=(omax_mmap>0)?omax_mmap<<1:2; // New (doubled) max_mmap
            // Create a new array twice the length of the old one
            ref_map* nmmap=new ref_map[nmax_mmap];
            // Copy the old array to the new one
            memcpy((void*)nmmap,(void*)fbase[k].mmap,sizeof(ref_map)*omax_mmap);
            nmmap[fbase[k].counter_mmap++]=ref_map(1000000.,1000000.,10000.,id,epsilon);
            delete [] fbase[k].mmap;
            fbase[k].mmap=nmmap;
            fbase[k].max_mmap=nmax_mmap;

            // Add to current list
            int mkid=match_base_id(k,id); // This would be different from the previous definition
            fbase[k].mmap[mkid].cc=ll;
            if(ll==1) {
                // First layer stores to the object
                of_list[id]->l.push_back(k);
                of_list[id]->l.push_back(id);
            } else {
                // Other layers
                l.push_back(k);
                l.push_back(id);
            }
            fbase[k].mmap[mkid].extrap_reset++;
        }
    }
    ov_buffer.clear();
}

inline void lbrmt_2d::set_layer(int i,int j,int k,int id) {
    // Extrapolation on one node
    double x=ax+i-1.5;
    double y=ay+j-1.5;
    int w=2;
    double nsolid;
    double sx;  // sum of x_i                
    double sy;  // sum of y_i
    double sxx; // sum of x_i^2
    double syy; // sum of y_i^2
    double sxy; // sum of x_i*y_i
    double sX;  // sum of X_i
    double sY;  // sum of Y_i
    double sxX; // sum of x_i*X_i
    double syX; // sum of y_i*X_i
    double sxY; // sum of x_i*Y_i
    double syY; // sum of y_i*Y_i
    double xavg;
    double yavg;
    double Xavg;
    double Yavg;    
    double ATA11;
    double ATA12;
    double ATA21;
    double ATA22;       
    double detATA;

    // Scan w*w window of solids
    // If nsolid<1, expand the window size
    int mkid=match_base_id(k,id);
    if(mkid==nobjs) {
        fprintf(stderr,"Can't match in set_layer(), t=%g, id=%d, i=%d, j=%d\n",time,id,i,j);
        exit(1);
    }
    int ll=fbase[k].mmap[mkid].cc;
    scan_solids(id,i,j,k,ll,x,y,w,nsolid,sx,sy,sxx,syy,sxy,sX,sY,sxX,syX,sxY,syY,xavg,yavg,Xavg,Yavg,ATA11,ATA12,ATA21,ATA22,detATA);

    // (A^T*A)^{-1} = [ATAI11 ATAI12, ATAI21 ATAI22]
    double ATAI11=1./detATA*ATA22;
    double ATAI12=-1./detATA*ATA12;
    double ATAI21=-1./detATA*ATA21;
    double ATAI22=1./detATA*ATA11;
    // sum of (x_i-xavg)*(X_i-Xavg)
    double xxXX=sxX-xavg*sX-Xavg*sx+nsolid*xavg*Xavg;
    // sum of (y_i-yavg)*(X_i-Xavg)
    double yyXX=syX-yavg*sX-Xavg*sy+nsolid*yavg*Xavg;
    // sum of (x_i-xavg)*(Y_i-Yavg)
    double xxYY=sxY-xavg*sY-Yavg*sx+nsolid*xavg*Yavg;
    // sum of (y_i-yavg)*(Y_i-Yavg)
    double yyYY=syY-yavg*sY-Yavg*sy+nsolid*yavg*Yavg;
    // Calculate x = [b, c]
    double bX=ATAI11*xxXX+ATAI12*yyXX;
    double cX=ATAI21*xxXX+ATAI22*yyXX;
    double bY=ATAI11*xxYY+ATAI12*yyYY;
    double cY=ATAI21*xxYY+ATAI22*yyYY;
    // Calculate extrapolated reference map value and update
    double extrapX=Xavg+bX*(x-xavg)+cX*(y-yavg);
    double extrapY=Yavg+bY*(x-xavg)+cY*(y-yavg);

    fbase[k].mmap[mkid].X=extrapX;
    fbase[k].mmap[mkid].Y=extrapY;
}

/** Scan the window for solids for reference map values extrapolation. */
inline void lbrmt_2d::scan_solids(int id,int i,int j,int k,int ll,double x,double y,int w,double &nsolid,double &sx,double &sy,
    double &sxx,double &syy,double &sxy,double &sX,double &sY,double &sxX,double &syX,double &sxY,double &syY,
    double &xavg,double &yavg,double &Xavg,double &Yavg,double &ATA11,double &ATA12,double &ATA21,double &ATA22,double &detATA) {
    // CHECK: prevent the code got stuck in this function
    int nsolve=0;
    int counter_solid=0;
    do {
        // CHECK: prevent code got stuck
        nsolve++;
        if(nsolve>10) printf("time = %g, nsolve = %d, i = %d, j = %d, ll = %d\n",time,nsolve,i,j,ll);
        // END OF CHECK
        w++; // Default w = 1+1 = 2
        nsolid=0.;
        sx=0.,sy=0.,sxx=0.,syy=0.,sxy=0.,sX=0.,sY=0.,sxX=0.,syX=0.,sxY=0.,syY=0.;
        // Clip the scan window to make sure it is not out of the simulation range
        int ilo=i-w<=0?0:i-w,ihi=i+w>=nxb?nxb-1:i+w,jlo=j-w<=0?0:j-w,jhi=j+w>=nyb?nyb-1:j+w;

        // Compute normal n_e
        double displacex=0.,displacey=0.; // Sum of displacement vectors of x_d - x_e
        for(int jj=jlo;jj<=jhi;jj++) {
            for(int ii=ilo;ii<=ihi;ii++) {
                // Find solid
                int kk=ii+jj*nxb;
                int mkkid=match_base_id(kk,id);    
                if(mkkid!=nobjs) {
                    int level=fbase[kk].mmap[mkkid].cc; // Level index in the extrapolation zone
                    if(level>=0&&level<ll) {
                        // x_e=(i,j), x_d=(ii,jj)
                        displacex+=ii-i;
                        displacey+=jj-j;
                    }
                }
            }
        }
        // Normal vector n_e
        double displacenorm=sqrt(displacex*displacex+displacey*displacey);
        double nex=-displacex/displacenorm;
        double ney=-displacey/displacenorm;

        // Weighted least-square extrapolation in the scan window
        for(int jj=jlo;jj<=jhi;jj++) {
            for(int ii=ilo;ii<=ihi;ii++) {
                // Find solid
                int kk=ii+jj*nxb;
                int mkkid=match_base_id(kk,id);    
                if(mkkid!=nobjs) {
                    int level=fbase[kk].mmap[mkkid].cc; // Level index in the extrapolation zone
                    if(level>=0&&level<ll) {
                        // Add coordinate-based weighting
                        // Weighting scheme in RMT3D to incorporate geometric information
                        static double powr[5]={1.,0.5,0.25,0.125,0.0625};
                        int idx=abs(jj-j)+abs(ii-i);
                        double r=idx<5?powr[idx]:pow(0.5,idx); // Change r when level<=2
                        if(level<=2) {
                            // Physical vector x'=x_e-x_d
                            // x_e=(i,j), x_d=(ii,jj)
                            double xpx=i-ii,xpy=j-jj;
                            double xpnorm=sqrt(xpx*xpx+xpy*xpy);
                            // Incorporate geometric when close to interface with n_e
                            double newr=(xpx*nex+xpy*ney)/xpnorm*r;
                            r=fmax(0.,newr);
                        }
                        nsolid+=r;
                        counter_solid+=1;
                        double X=fbase[kk].mmap[mkkid].X,Y=fbase[kk].mmap[mkkid].Y;
                        double xx=x+ii-i;
                        double yy=y+jj-j;
                        sx+=xx*r;
                        sy+=yy*r;
                        sxx+=xx*xx*r;
                        syy+=yy*yy*r;
                        sxy+=xx*yy*r;
                        sX+=X*r;
                        sY+=Y*r;
                        sxX+=xx*X*r;
                        syX+=yy*X*r;
                        sxY+=xx*Y*r;
                        syY+=yy*Y*r;
                    }
                }
            }
        }
        xavg=sx/nsolid;
        yavg=sy/nsolid;
        Xavg=sX/nsolid;
        Yavg=sY/nsolid;     
        // A^T*A = [ATA11 ATA12, ATA21 ATA22]
        ATA11=sxx-sx*sx/nsolid;
        ATA12=sxy-sx*sy/nsolid;
        ATA21=ATA12;
        ATA22=syy-sy*sy/nsolid;
        detATA=ATA11*ATA22-ATA12*ATA21;
    } while (nsolid<1. && counter_solid<3 && (fabs(detATA)<1e-12));
}

/** Calculate the deviatoric part of an incompressible solid stress tensor.
 * \param[in] (Xx,Xy,Yx,Yy) the components of the Jacobian.
 * \param[out] (J,BP11,BP12,BP21,BP22) the determinant of F and the components of deviatoric B'. */
inline void lbrmt_2d::calc_fft(double X,double Y,double Xx,double Xy,double Yx,double Yy,double &J,double &BP11,double &BP12,double &BP21,double &BP22,object *o) {
    // Calculate the deformation gradient
    // F = [F11 F12, F21 F22]
    double det=Xx*Yy-Xy*Yx;
    double F11=1./det*Yy;
    double F12=-1./det*Xy;
    double F21=-1./det*Yx;
    double F22=1./det*Xx;
    // Apply the active deformation F_a on the deformation gradient tensor
    double xx,yy;
    o->transform(X,Y,xx,yy); // (rmap.X,rmap.Y) world space, (xx,yy) object space
    o->actuate(xx,yy,F11,F12,F21,F22);
    // Calculate the determinant of F
    // J = det(F)
    J=F11*F22-F12*F21;
    // Calculate the transpose of the deformation gradient
    // FT = [FT11 FT12, FT21 FT22]
    double FT11=F11;
    double FT12=F21;
    double FT21=F12;
    double FT22=F22;
    // Calculate FFT
    // FFT = F * F^T = B
    double FFT11=F11*FT11+F12*FT21;
    double FFT12=F11*FT12+F12*FT22;
    double FFT21=F21*FT11+F22*FT21;
    double FFT22=F21*FT12+F22*FT22;
    double traceFFT=FFT11+FFT22;
    // Compute the deviatoric part of a tensor
    // B' = B - 1/3 * (trace(B)+1) * I
    // BP = [BP11 BP12, BP21 BP22]
    BP11=FFT11-1./3*(traceFFT+1.0);
    BP12=FFT12;
    BP21=FFT21;
    BP22=FFT22-1./3*(traceFFT+1.0);
}

/** Compute the determinant of the deformation gradient.
 * \param[in] (Xx,Xy,Yx,Yy) the components of the Jacobian.
 * \param[out] J the determinant of F. */
inline double lbrmt_2d::calc_detF(double Xx,double Xy,double Yx,double Yy) {
    // Calculate the deformation gradient
    // F = [F11 F12, F21 F22]
    double det=Xx*Yy-Xy*Yx;
    double F11=1./det*Yy;
    double F12=-1./det*Xy;
    double F21=-1./det*Yx;
    double F22=1./det*Xx;
    // Calculate the determinant of F
    // J = det(F)
    double J=F11*F22-F12*F21;
    return J;
}

/** Compute the inverse of determinant of the deformation gradient.
 * \param[in] (Xx,Xy,Yx,Yy) the components of the Jacobian.
 * \param[out] det the determinant of F. */
inline double lbrmt_2d::calc_invdetF(double Xx,double Xy,double Yx,double Yy) {
    // Calculate the deformation gradient
    // F = [F11 F12, F21 F22]
    double det=Xx*Yy-Xy*Yx;
    return det;
}

/** Compute contact force.
 * \param[in] phi the level set value.
 * \return the contact force. */
inline double lbrmt_2d::contact_force(double phi) {
    return phi>=epsilon?0.:0.5*(1.-phi/epsilon);
}

/** Return the ref_map object index in the mmap that matches the object ID.
 * \param[in] i the x index.
 * \param[in] j the y index.
 * \param[in] k the index of the f array.
 * \param[in] id the object ID.
 * \return the ref_map object index in the mmap. */
inline int lbrmt_2d::match_id(int i,int j,int k,int id) {
    for(int q=0;q<f[k].counter_mmap;q++) {
        if(f[k].mmap[q].id==id) return q;
    }
    return nobjs;
}

/** Return the ref_map object index in the mmap that matches the object ID for fbase.
 * \param[in] k the index of the fbase array.
 * \param[in] id the object ID.
 * \return the ref_map object index in the mmap. */
inline int lbrmt_2d::match_base_id(int k,int id) {
    for(int q=0;q<fbase[k].counter_mmap;q++) {
        if(fbase[k].mmap[q].id==id) return q;
    }
    // TODO: better way to detect not matching
    return nobjs;
}

/** A smoothed Heaviside function with a transition region.
 * \param[in] phi the level set value.
 * \return the smoothed Heaviside value. */
inline double lbrmt_2d::heaviside(double phi) {
    if(phi>=epsilon) return 0.;       // Fluid
    else if(phi<=-epsilon) return 1.; // Solid
    else {
        // 2020 JFM Heaviside function
        double xx=phi/epsilon;
        return 1.-0.5*(1.+xx+1./M_PI*sin(M_PI*xx));
    }
}

/** A smoothed delta function with a transition region.
 * \param[in] phi the level set value.
 * \return the smoothed delta value. */
inline double lbrmt_2d::delta(double phi) {
    if(phi>epsilon) return 0.;       // Fluid
    else if(phi<-epsilon) return 0.; // Solid
    else {
        // 2020 JFM delta function
        double xx=phi/epsilon;
        return 0.5*(1.+cos(M_PI*xx));
    }
}

/** Count total number of ref_map objects. */
void lbrmt_2d::diag_count_ref_map() {
    int count_ref_map=0;
    int count_max_mmap=0;
#pragma omp parallel for
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;
            // Check max mmap length
            int f_max_mmap=f[k].max_mmap;
            count_max_mmap+=f_max_mmap;
            // Check number of ref_map objects
            int f_counter_mmap=f[k].counter_mmap;
            if(f_counter_mmap>0)
            count_ref_map+=f_counter_mmap;
            
        }
    }
    if(int(time)%1000==0) printf("time = %g, total ref_map = %d, max mmap = %d\n",time,count_ref_map,count_max_mmap);
}

/** Compute the total density of the simulation domain. */
void lbrmt_2d::diag_total_rho() {
    double total_rho=0.;
#pragma omp parallel for reduction(+:total_rho)
    for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
            int k=j*nxb+i;
            total_rho+=f[k].rho;
        }
    }
    if(int(time)%100==0) printf("time = %g, total density = %.16f, average density = %.16f\n",time,total_rho,total_rho/nx/ny);
}
