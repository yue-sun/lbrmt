#include <cstdio>
#include <cstring>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include "convert_cfg.hh"
#include "convert_unit.hh"
#include "lbrmt_2d.hh"
#include "omp.h"

int main(int argc,char **argv) {
    // Check for the correct number of command-line arguments
    if(argc!=2&&(argc!=3||strcmp(argv[1],"-m")!=0)) {
        fputs("Syntax: ./sim_fsi_mix [-m] <input_file>\n\n"
              "-m pass Mach number limit\n"
              "The input file should have a '.cfg' suffix\n",stderr);
        return 1;
    }

    // Read in the parameters from the input file
    char* infile=argv[argc-1];
    char* outdir=NULL;       // Output directory name
    char* outdir_mmap=NULL;  // Output multimaps directory name
    int niters;              // Number of total simulation iterations
    int nout;                // Number of output frames
    int nx;                  // Number of grid points in x direction
    int ny;                  // Number of grid points in y direction
    double epsilon;          // Blur zone half-width
    int bctype;              // Fluid boundary conditions
    double rho_f;            // Dimensionless fluid density
    double Re;               // Reynolds number
    double tau;              // LB relaxation time
    double nu_ph;            // Physical fluid kinematic viscosity
    double ux0;              // Initial fluid velocity in x direction
    double uy0;              // Initial fluid velocity in y direction
    double fx_fluid_ph;      // Physical fluid force density in x direction
    double fy_fluid_ph;      // Physical fluid force density in y direction
    double dx;               // Physical grid spacing
    double dt;               // Physical timestep
    double C_rho;            // Density scale (benchmark dimensionless to physical)
    int nobjs;               // Number of solid objects
    convert_cfg(infile,outdir,outdir_mmap,niters,nout,nx,ny,epsilon,bctype,
        rho_f,Re,tau,nu_ph,ux0,uy0,fx_fluid_ph,fy_fluid_ph,dx,dt,C_rho,nobjs);

    // Set up LB simulation parameters
    double L=fmin(nx,ny);           // Characteristic length scale
    double c=sqrt(1./3);            // LB speed of sound
    double cs2=1./3;                // LB speed of sound squared
    double nu=cs2*(tau-0.5);        // LB kinematic viscosity

    // Determine dt or tau based on input
    // nu_ph = cs^2 * (tau - 0.5) dx^2 / dt
    if(tau==0.) {
        // Determine tau based on dt and nu_ph
        nu=nu_ph/dx/dx*dt;        // LB kinematic viscosity
        tau=nu/cs2+0.5;
    } else if (dt==0. && tau>0.5) {
        // Determine dt based on tau and nu_ph
        dt=nu/nu_ph*dx*dx;        // Physical timestep
    }
    
    double umax_ph=Re*nu_ph/(L*dx); // Physical umax in fluid
    double umax=Re*nu/L;            // LB umax in fluid
    double Ma=umax/c;               // LB Mach number

    if(Ma>0.3 && strcmp(argv[1],"-m")!=0) {
        printf("Error: Compressible Mach number %g is too large (Ma*>0.3).\n",Ma);
        return -1;
    }
    if(tau<=0.5 && strcmp(argv[1],"-m")!=0) {
        printf("Error: Relaxation time tau = %g is too small (tau*>0.5).\n",tau);
        return -1;
    }

    // Convert parameter units and calculate LBM values
    double fx_fluid=0.; // LB fluid force density in x direction
    double fy_fluid=0.; // LB fluid force density in y direction
    if(bctype==0) {
        // Arbitrary external force for Poiseuille flow
        double C_g=dx/dt/dt; // Conversion factor of force density
        fx_fluid+=convert_force(fx_fluid_ph,C_g);
        // Poiseuille flow with the given Reynolds number (i.e. umax)
        fx_fluid+=umax*8.*nu/L/L;
    }

    // Get bctype text string
    std::string bctext;
    switch(bctype) {
        case 0: bctext="periodic channel";break;
        case 3: bctext="fully no-slip box";break;
        case 4: bctext="fully periodic box";break;
        case 7: bctext="lid-driven cavity";break;
        case 8: bctext="Taylor–Green vortex";break;
    }

    // Print simulation parameters
    printf("Input file: %s\n"
           "Output directory: %s\n"
           "==========================================================================\n"
           "Parameters                              Variables      Values\n"
           "==========================================================================\n"
           "total simulation iterations             niters:        %d\n"
           "number of output frames                 nout:          %d\n"
           "number of grid points in x              nx:            %d\n"
           "number of grid points in y              ny:            %d\n"
           "blur zone half-width                    epsilon:       %g\n"
           "fluid boundary conditions               bctype:        %s\n"
           "dimensionless fluid density             rho_f:         %g\n"
           "initial fluid velocity in x             ux0:           %g\n"
           "initial fluid velocity in y             uy0:           %g\n"
           "==========================================================================\n"
           "Reynolds number                         Re:            %g\n"
           "LB relaxation time                      tau:           %g\n"
           "physical fluid kinematic viscosity      nu_ph:         %g\n"
           "LB fluid kinematic viscosity            nu:            %g\n"
           "physical grid spacing                   dx:            %g\n"
           "physical time spacing in fluid          dt:            %g\n"
           "physical umax in fluid                  umax_ph:       %g\n"
           "LB umax in fluid                        umax:          %g\n"
           "LB Mach number                          Ma:            %g\n"
           "LB fluid force density in x             fx_fluid:      %g\n"
           "LB fluid force density in y             fy_fluid:      %g\n"
           "==========================================================================\n"
           "number of solid objects                 nobjs:         %d\n"
           "==========================================================================\n",
           infile,outdir,niters,nout,nx,ny,epsilon,bctext.c_str(),
           rho_f,ux0,uy0,Re,tau,nu_ph,nu,dx,dt,umax_ph,umax,Ma,fx_fluid,fy_fluid,
           nobjs);
    
    // Set up solid geometry and physical properties
    std::vector<object*> objs;
    std::vector<mat_const*> mc;
    if(nobjs>0) {
        // Read in the parameters from the input file for solid objects
        int obj_type; // Solid geometry type
        double rho_s; // Dimensionless solid density
        double G;     // Dimensionless solid shear modulus
        double cx;    // Dimensionless center of solid object in x direction
        double cy;    // Dimensionless center of solid object in y direction
        double cr;    // Dimensionless size of solid object
        convert_solid_cfg(infile,obj_type,rho_s,G,cx,cy,cr);

        double C_G=C_rho*dx*dx/dt/dt; // Conversion factor of small-strain shear modulus
        bool set_velocity=false;

        // Check timestep stability
        double dt1=sqrtl(rho_s/G)*dx;
        if(dt>dt1 && strcmp(argv[1],"-m")!=0) {
            printf("Error: Timestep is larger than shear wave limit (dt_I = %g).\n",dt1);
            return -1;
        }
        
        object *os;
        mat_const *one_mc;
        // Create object based on solid type and get obj_type text string
        std::string objtext;

        // Four rotors
        if(nobjs==4) {
            for(int k=-1;k<1;k++) {
                for(int l=-1;l<1;l++) {
                    switch(obj_type) {
                        case 5: os=new obj_smooth_rotor(cx*nx+0.4*nx*(l+0.5),cy*ny+0.4*nx*(k+0.5),cr*nx,4.*cr*nx,cr*nx,0.5*cr*nx,k%2==0?0.:M_PI/3.,k%2==0?M_PI:-M_PI,2*M_PI/(niters));
                            objtext="smooth rotor";
                            break;
                    }
                    objs.push_back(os);
                    // Assign material constants
                    one_mc=new mat_const(convert_G(G,C_G),rho_s,set_velocity);
                    mc.push_back(one_mc);
                }
            }
        }

        // Print simulation parameters for solid objects
        printf("shear wave timestep limit               dt_I:          %g\n"
               "solid geometry type                     obj_type:      %s\n"
               "dimensionless solid density             rho_s          %g\n"
               "dimensionless solid shear modulus       G:             %g\n"
               "length scale of solid object            cr:            %g\n"
               "==========================================================================\n",
               dt1,objtext.c_str(),rho_s,G,cr);
    }

    // Create the solver
    lbrmt_2d fsi(rho_f,Re,tau,epsilon,nx,ny,0.,nx,0.,ny,bctype,nobjs,outdir,outdir_mmap);
    // Add solid objects
    for(int k=0;k<nobjs;k++) fsi.add_object(objs[k],mc[k],k);
    // Initialize
    fsi.initialize(rho_f,ux0,uy0,bctype,fx_fluid,fy_fluid);
    // Start measuring time
    printf("Start simulation:\n");
    struct timespec wall_begin,wall_end;
    clock_gettime(CLOCK_REALTIME, &wall_begin);                       // Begin wall time
    clock_t cpu_begin=clock();                                        // Begin CPU clock tick
    // Run simulaton
    fsi.solve(niters,nout+1,bctype,fx_fluid,fy_fluid);
    // Run simulaton with timing
    // fsi.solve_time(niters,nout+1,bctype,fx_fluid,fy_fluid);
    // Stop measuring time and calculate the elapsed time
    clock_gettime(CLOCK_REALTIME,&wall_end);                          // End wall time
    clock_t cpu_end=clock();                                          // End CPU clock tick
    long seconds=wall_end.tv_sec-wall_begin.tv_sec;
    long nanoseconds=wall_end.tv_nsec-wall_begin.tv_nsec;
    double wall_time=seconds+nanoseconds*1e-9;                        // Wall time
    double cpu_time=double(cpu_end-cpu_begin)/double(CLOCKS_PER_SEC); // CPU time
    int num_threads;
#pragma omp parallel
    {
        num_threads=omp_get_num_threads();
    }
           
    printf("End simulation successfully!\n"
           "==========================================================================\n"
           "Number of threads: %d\n"
           "Wall time: %g\n"
           "CPU time: %g\n",
           num_threads,wall_time,cpu_time);

    /** Free dynamically allocated objects. */
    for(int k=0;k<nobjs;k++) {
        delete objs[k];
        delete mc[k];
    }
    delete [] outdir;
    delete [] outdir_mmap;
}