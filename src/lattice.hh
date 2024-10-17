#ifndef LATTICE_HH
#define LATTICE_HH

#include <cmath>
#include <cstdio>
#include <vector>

#include "ref_map.hh"

// LB sound speed constants
static const double cs2=1./3,cs4=1./9;
// Weights for D2Q9 model
static const double w0=4./9,w1=1./9,w2=1./36;

class lattice {
    public:
        /***************** Lattice Boltzmann method variables. *****************/
        /** Macroscopic hydro variables. */
        double rho,ux,uy;
        /** Populations at current state. */
        double f0,f1,f2,f3,f4,f5,f6,f7,f8;
        /** Equilibrium populations (also used to store streamed populations). */
        double feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8;
        /** Force density at current state. */
        double F0,F1,F2,F3,F4,F5,F6,F7,F8;
        /** Smooth flux correction on populations at current state. */
        double fc0,fc1,fc2,fc3,fc4,fc5,fc6,fc7,fc8;
        /** Net (smooth Heavisided solid and fluid) external force density. */
        double Fx,Fy;
        /** Magnitude of solid stress. */
        double ss;
        /** Obstacle flag. */
        bool obs;

        /***************** Reference map technique variables. ******************/
        /** A multimaps pointer to all ref_map objects. */
        ref_map* mmap;
        /** Counter of the mmap array, i.e. number of current ref_map objects. */
        int counter_mmap;
        /** Max length of the mmap array, can be continuously doubled. */
        int max_mmap;
        /** Components of solid stress. */
        double sla,slb,slc,sld; // solid stress left
        double sba,sbb,sbc,sbd; // solid stress right
        /** Total solid fraction. */
        double total_lambda;
        /** Left and bottom determinant of inverse of deformation gradient. */
        double invdetF_left,invdetF_bottom;
        /** Total determinant of deformation gradient. */
        double total_detF;
        /** Average of the inverse of detF. */
        double avg_invdetF;

        /** Destroy function. */
        virtual ~lattice() {}

        /** Initialize hydro variables with force fields, p.p 249 Kruger's book.
         * \param[in] macro_rho the Heavisided macroscopic density.
         * \param[in] macro_ux the x component of macroscopic velocity.
         * \param[in] macro_uy the y component of macroscopic velocity.
         * \param[in] fx_h the x component of Heavisided external force density.
         * \param[in] fy_h the y component of Heavisided external force density. */
        inline void init_hydro(double macro_rho,double macro_ux,double macro_uy,double fx_h,double fy_h) {
            rho=macro_rho;
            ux=macro_ux-0.5/rho*fx_h;
            uy=macro_uy-0.5/rho*fy_h;
        }

        /** Initialize hydro variables to equilibrium. */
        inline void init_pop() {
            f0=feq0;
            f1=feq1;
            f2=feq2;
            f3=feq3;
            f4=feq4;
            f5=feq5;
            f6=feq6;
            f7=feq7;
            f8=feq8;
        }

        /** Calculate updated macroscopic hydro variables rho, ux and uy with the 
         * smoothed force densities, p.p 233 Kruger's book. */
        inline void hydro() {
            rho=f0+f1+f2+f3+f4+f5+f6+f7+f8;
            ux=1./rho*(f1+f5+f8-f3-f6-f7)+0.5/rho*Fx;
            uy=1./rho*(f2+f5+f6-f4-f7-f8)+0.5/rho*Fy;
        }

        /** Calculate the equilibrium populations.
         * \param[in] drho target density difference. */
        inline void equilibrium(double drho) {
            double uxsq=ux*ux;
            double uysq=uy*uy;
            double usq=uxsq+uysq;

            // Smooth flux correction
            fc1=w1*drho;
            fc2=w1*drho;
            fc3=w1*drho;
            fc4=w1*drho;
            fc5=w2*drho;
            fc6=w2*drho;
            fc7=w2*drho;
            fc8=w2*drho;
            fc0=-fc1-fc2-fc3-fc4-fc5-fc6-fc7-fc8;

            feq0=w0*rho*(1.-0.5*usq/cs2)                                       ;
            feq1=w1*rho*(1.+ux/cs2+0.5*uxsq/cs4-0.5*usq/cs2)                   ;
            feq2=w1*rho*(1.+uy/cs2+0.5*uysq/cs4-0.5*usq/cs2)                   ;
            feq3=w1*rho*(1.-ux/cs2+0.5*uxsq/cs4-0.5*usq/cs2)                   ;
            feq4=w1*rho*(1.-uy/cs2+0.5*uysq/cs4-0.5*usq/cs2)                   ;
            feq5=w2*rho*(1.+( ux+uy)/cs2+0.5*( ux+uy)*( ux+uy)/cs4-0.5*usq/cs2);
            feq6=w2*rho*(1.+(-ux+uy)/cs2+0.5*(-ux+uy)*(-ux+uy)/cs4-0.5*usq/cs2);
            feq7=w2*rho*(1.+(-ux-uy)/cs2+0.5*(-ux-uy)*(-ux-uy)/cs4-0.5*usq/cs2);
            feq8=w2*rho*(1.+( ux-uy)/cs2+0.5*( ux-uy)*( ux-uy)/cs4-0.5*usq/cs2);
        }

        /** Collide to calculate the newly-updated population for the next timestep.
         * \param[in] invtau 1/tau. */
        inline void collide(double invtau) {
            // Collision operator
            double Omega0=-invtau*(f0-feq0)-invtau*fc0;
            double Omega1=-invtau*(f1-feq1)-invtau*fc1;
            double Omega2=-invtau*(f2-feq2)-invtau*fc2;
            double Omega3=-invtau*(f3-feq3)-invtau*fc3;
            double Omega4=-invtau*(f4-feq4)-invtau*fc4;
            double Omega5=-invtau*(f5-feq5)-invtau*fc5;
            double Omega6=-invtau*(f6-feq6)-invtau*fc6;
            double Omega7=-invtau*(f7-feq7)-invtau*fc7;
            double Omega8=-invtau*(f8-feq8)-invtau*fc8;

            // Update populations with collision operator
            f0+=Omega0;
            f1+=Omega1;
            f2+=Omega2;
            f3+=Omega3;
            f4+=Omega4;
            f5+=Omega5;
            f6+=Omega6;
            f7+=Omega7;
            f8+=Omega8;

            // Update populations with source term
            double inv2tau=0.5*invtau;
            f0+=(1.-inv2tau)*F0;
            f1+=(1.-inv2tau)*F1;
            f2+=(1.-inv2tau)*F2;
            f3+=(1.-inv2tau)*F3;
            f4+=(1.-inv2tau)*F4;
            f5+=(1.-inv2tau)*F5;
            f6+=(1.-inv2tau)*F6;
            f7+=(1.-inv2tau)*F7;
            f8+=(1.-inv2tau)*F8;
        }

        /** Add the net external force density with Guo's forcing scheme, p.p 236-237 Kruger's book.
         * \param[in] tau the relaxation time. */
        inline void force(double tau) {
            double Fxux=Fx*ux;
            double Fyuy=Fy*uy;
            double Fdotu=Fxux+Fyuy;

            // Second-order force discretization in velocity space
            F0=w0*(-Fdotu/cs2);
            F1=w1*( Fx/cs2+Fxux/cs4-Fdotu/cs2);
            F2=w1*( Fy/cs2+Fyuy/cs4-Fdotu/cs2);
            F3=w1*(-Fx/cs2+Fxux/cs4-Fdotu/cs2);
            F4=w1*(-Fy/cs2+Fyuy/cs4-Fdotu/cs2);
            F5=w2*(( Fx+Fy)/cs2+( Fx+Fy)*( ux+uy)/cs4-Fdotu/cs2);
            F6=w2*((-Fx+Fy)/cs2+(-Fx+Fy)*(-ux+uy)/cs4-Fdotu/cs2);
            F7=w2*((-Fx-Fy)/cs2+(-Fx-Fy)*(-ux-uy)/cs4-Fdotu/cs2);
            F8=w2*(( Fx-Fy)/cs2+( Fx-Fy)*( ux-uy)/cs4-Fdotu/cs2);
        }

        /** Update new populations in streaming. */
        inline void update_f() {
            init_pop();
        }
};

#endif