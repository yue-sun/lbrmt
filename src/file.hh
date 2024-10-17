/** This C++ code is from Chris Rycroft AM225 sample code to organize file I/O better.
 * https://github.com/chr1shr/am225_examples/blob/master/5a_fluid_sim/common.cc
 * https://github.com/chr1shr/am225_examples/blob/master/5a_fluid_sim/common.hh */

#ifndef FILE_HH
#define FILE_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>

/** Open a file, and check the return value to ensure that the operation was successful. */
FILE* safe_fopen(const char* filename,const char* mode);

/** Read data from a file and check the return value to ensure that the operation was successful.
 * If not successful, print an error message and exit. */
void fatal_error(const char *p,int code);

/** Function for printing fatal error messages and exiting. */
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p);

#endif