/** This C++ code is adapted from Chris Rycroft VDropImpact sample code to
 * parse simulation config file into the solver.
 * https://github.com/chr1shr/vdropimpact/blob/master/fileinfo.cc
 * https://github.com/chr1shr/vdropimpact/blob/master/fileinfo.hh */

#ifndef CONVERT_CFG_HH
#define CONVERT_CFG_HH

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

/** The size of the temporary buffer for parsing the input file. */
const int cfg_buf_size=512;

/** Load .cfg file and convert input texts into simulation parameters. */
void convert_cfg(const char* infile,char* &outdir,char* &outdir_mmap,
    int &niters,int &nout,int &nx,int &ny,double &epsilon,int &bctype,
    double &rho_f,double &Re,double &tau,double &nu_ph,
    double &ux0,double &uy0,double &fx_fluid_ph,double &fy_fluid_ph,
    double &dx, double &dt,double &C_rho,int &nobjs);

/** Load .cfg file and convert input texts into simulation parameters.
 * for solid objects. */
void convert_solid_cfg(const char* infile,int &obj_type,double &rho_s,
    double &G,double &cx,double &cy,double &cr);

/** Find the next token in a string and interpret it as a double precision
 * floating point number. If none is availble, it gives an error message. */
inline double next_double(int ln);

/** Find the next token in a string, interpret it as a double precision 
 * floating point number, and checks that there are no subsequent values. */
inline double final_double(int ln);

/** Find the next token in a string, interpret it as an integer, and check
 * that there are no subsequent values. */
inline int final_int(int ln);

/** Test to see if two strings are equal. */
inline bool se(const char *p1,const char *p2);

/** Find the next token in a string, and if none is availble, raise error. */
char* next_token(int ln);

/** Check that there are no subsequent values. */
void check_no_more(int ln);

/** Check that a parameter is a valid positive value. */
void check_invalid(double val,const char *p);

/** Generate a random double between two doubles. */
double rand_ab(double a,double b);

#endif