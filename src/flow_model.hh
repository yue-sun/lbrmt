#ifndef FLOW_MODEL
#define FLOW_MODEL

#include <cmath>

#include "lbrmt_2d.hh"

/** Compute the exact solution for Taylor–Green vortex decay.
 * \param[in] t simulation time.
 * \param[in] (x,y) the current position.
 * \param[in] scale scaling factor for convergence study.
 * \param[in,out] rho fluid density.
 * \param[in,out] (ux,uy) the velocity. */
void lbrmt_2d::taylor_green(double t,double x,double y,int scale,double &rho,double &ux,double &uy) {
    // Set constants
    double u0=0.04/scale;               // Initial velocity scale
    double kx=2.*M_PI/nx,ky=2.*M_PI/ny; // Components of wave factor
    double td=1./nu/(kx*kx+ky*ky);      // Vortex decay time

    // Compute velocity
    ux=-u0*sqrt(ky/kx)*cos(kx*x)*sin(ky*y)*exp(-1.*t/td);
    uy= u0*sqrt(kx/ky)*sin(kx*x)*cos(ky*y)*exp(-1.*t/td);
    
    // Compute pressure and density
    double p=-0.25*rho_f*u0*u0*(ky/kx*cos(2.*kx*x)+kx/ky*cos(2.*ky*y))*exp(-2.*t/td);
    rho=rho_f+3.*p;
}

#endif