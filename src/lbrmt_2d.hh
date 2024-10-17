#ifndef LBRMT_2D_HH
#define LBRMT_2D_HH

#include <vector>

#include "file.hh"
#include "lattice.hh"
#include "obj_field.hh"
#include "object.hh"

class lbrmt_2d {
    public:
        /******************************* Physical parameters *********************************/
        /** Fluid density. */
        const double rho_f;
        /** Reynolds number. */
        const double Re;
        /** Relaxation time. */
        const double tau;
        /** Kinematic viscosity. */
        const double nu;
        /** Blur zone half-width. */
        const double epsilon;

        /*************************** Simulation domain dimensions ****************************/
        /** Number of grid points in the x direction, and that padded with buffer nodes. */
        const int nx,nxb;
        /** Number of grid points in the y direction, and that padded with buffer nodes. */
        const int ny,nyb;
        /** Lower and upper bound of simulation domain in the x direction. */
        const double ax,bx;
        /** Lower and upper bound of simulation domain in the y direction. */
        const double ay,by;

        /*********************************** Miscellaneous ***********************************/
        /** Boundary condition types. */
        int bctype;
        /** Initial fluid velocity. */
        double ux0,uy0;
        /** Number of all solid objects. */
        int nobjs;
        /** Simulation time. */
        double time;
        /** Overflow buffer for un-instantiated ref_map objects ID and (i,j). */
        std::vector<int> ov_buffer;
        /** Counter of overflow buffer. */
        int counter_ov;

        /********************************* Simulation arrays *********************************/
        /** A pointer to all nodes including the buffer nodes. */
        lattice *fbase;
        /** A pointer to only physically-defined nodes. */
        lattice *f;
        /** A pointer to a boolean array of rigid obstacles. */
        bool *obs;

        /** Instantitate class object. */
        lbrmt_2d(double rho_f_,double Re_,double tau_,double epsilon_,
            int nx_,int ny_,double ax_,double bx_,double ay_,double by_,
            int bctype_,int nobjs, char *outdir_, char *outdir_mmap_,
            const char *obsfile_=NULL,double ux0=0.,double uy0=0.);
        /** Destroy class object. */
        ~lbrmt_2d();

        /****************************** Initialization routine *******************************/
        /** Add solid objects. */
        void add_object(object *o,mat_const *mc,int oid);
        /** Initialize reference map of solid nodes. */
        void init_remap();
        /** Initialize hydro variables of all nodes. */
        void init_hydro(double rho,double macro_ux,double macro_uy,double fx_fluid=0.,double fy_fluid=0.);
        /** Initialize populations of all nodes. */
        void init_pop();
        /** Initialization function called in the main routine. */
        void initialize(double rho_f,double macro_ux,double macro_uy,int bctype,double fx_fluid=0.,double fy_fluid=0.);

        /*********************************** Main routine ************************************/
        /** Calculate updated hydro variables. */
        void hydro();
        /** Calculate updated reference map field. */
        void remap();
        /** Extrapolate the reference map values in the blur zone positive direction. */
        void extrap();
        /** Calculate the total solid fractions for all objects. */
        void total_solid_fraction();
        /** Calculate all stress for all solid nodes. */
        void compute_stress();
        /** Helper function to calculate penalized-incompressible solid stress tensor. */
        template<bool left>
        void solid_stress(int i,int j,int k,int id,ref_map &rmap_center,object *o,double G);
        /** Helper function to calculate collision stress based on contact forces between solid objects. */
        template<bool left>
        void collision_stress(int i,int j,int k,int id,int q,ref_map &rmap_center,object *o,double G);
        /** Calculate the equilibrium populations. */
        void equilibrium();
        /** Collide and calculate the populations for the next timestep. */
        void collide();
        /** Compute the external force density. */
        void force_density(double fx_fluid,double fy_fluid,int bctype);
        /** Add forcing term. */
        void force();
        /** Set the boundary conditions. */
        void bc(int bctype);
        /** Move populations along D2Q9 directions. */
        void stream();
        /** Add obstacle. */
        void obstacle();
        /** Main solve routine. */
        void solve(int niters,int nout,int bctype,double fx_fluid=0.,double fy_fluid=0.);
        /** Main solve routine with timing. */
        void solve_time(int niters,int nout,int bctype,double fx_fluid=0.,double fy_fluid=0.);

        /**************************** Boundary conditions routine ****************************/
        /** Periodic channel boundary conditions for lbm variables. */
        void pbc_lbm();
        /** Open channel boundary conditions for lbm variables. */
        void obc_lbm();
        /** Fully no-slip box boundary conditions for lbm variables. */
        void fnbc_lbm();
        /** Fully periodic channel boundary conditions for lbm variables. */
        void fpbc_lbm();
        /** Lid-driven cavity boundary conditions for lbm variables. */
        void ldc_lbm();
        /** Taylor–Green vortex for lbm variables. */
        void taylor_green(double t,double x,double y,int scale,double &rho,double &ux,double &uy);
        
        /***************************** Post-processing routine *******************************/
        /** Set output directory. */
        void set_dir(char *outdir_,char *outdir_mmap_);
        /** Write to file. */
        void output(int fr);

        /******************************** Diagnostic routine *********************************/
        /** Count total number of ref_map objects. */
        void diag_count_ref_map();
        /** Compute the total density of the simulation domain. */
        void diag_total_rho();

    private:
        /** An array of pointers to obj_field classes, which contain the fields for each solid
         *  in the simulation. */
        obj_field **of_list;
        /** Output directory filename. */
        char *outdir;
        /** Output directory filename for multimaps. */
        char *outdir_mmap;
        /** Buffer for assembling output filename. */
        char *outbuf;
        /** Buffer for assembling output filename for multimaps. */
        char *outbuf_mmap;
        /** Pointer to output file. */
        FILE *fp;
        /** Pointer to output file for multimaps. */
        FILE *fp_mmap;
        /** Second-order ENO upwinding method. */
        inline double eno2(double v0,double v1,double v2,double v3);
        /** Build each layer in the extrapolation zone. */
        inline void build_layer(int i,int j,int k,int ll,std::vector<int> &l,int id);
        /** Set the extrapolated reference map value of each layer in the blur zone positive direction. */
        inline void set_layer(int i,int j,int k,int id);
        /** Scan the window for solids for reference map values extrapolation. */
        inline void scan_solids(int id,int i,int j,int k,int ll,double x,double y,int w,double &nsolid,double &sx,double &sy,double &sxx,double &syy,double &sxy,
            double &sX,double &sY,double &sxX,double &syX,double &sxY,double &syY,double &xavg,double &yavg,double &Xavg,double &Yavg,
            double &ATA11,double &ATA12,double &ATA21,double &ATA22,double &detATA);
        /** Resolve overflow buffer in extrapolation. */
        inline void resolve_overflow(int ll,std::vector<int> &l);
        /** Calculate the deviatoric part of an incompressible solid stress tensor. */
        inline void calc_fft(double X,double Y,double Xx,double Xy,double Yx,double Yy,double &J,double &BP11,double &BP12,double &BP21,double &BP22,object *o);
        /** Calculate the determinant of the deformation gradient tensor. */
        inline double calc_detF(double Xx,double Xy,double Yx,double Yy);
        /** Calculate the inverse of determinant of the deformation gradient tensor. */
        inline double calc_invdetF(double Xx,double Xy,double Yx,double Yy);
        /** Compute contact force. */
        inline double contact_force(double phi);
        /** Return the ref_map object index in the mmap that matches the object ID. */
        inline int match_id(int i,int j,int k,int id);
        /** Return the ref_map object index in the mmap that matches the object ID for fbase. */
        inline int match_base_id(int k,int id);
        /** A smoothed Heaviside function with a transition region. */
        inline double heaviside(double phi);
        /** A smoothed delta function with a transition region. */
        inline double delta(double phi);
        /** Wall boundary force. */
        inline double delta_func_wall(const double& x) {
            return (1.-x)/x;
        }
};

#endif