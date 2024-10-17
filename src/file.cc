/** This C++ code is from Chris Rycroft AM225 sample code to organize file I/O better.
 * https://github.com/chr1shr/am225_examples/blob/master/5a_fluid_sim/common.cc
 * https://github.com/chr1shr/am225_examples/blob/master/5a_fluid_sim/common.hh */

#include "file.hh"

/** Open a file, and check the return value to ensure that the operation was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return the file handle. */
FILE* safe_fopen(const char* filename,const char* mode) {
    FILE *temp=fopen(filename,mode);
    if(temp==NULL) fprintf(stderr,"Error opening file \"%s\"\n",filename);
    return temp;
}

/** Read data from a file and check the return value to ensure that the operation was successful.
 * If not successful, print an error message and exit.
 * \param[in] ptr the memory to write to.
 * \param[in] size the size of each element to be read.
 * \param[in] count the number of elements to be read.
 * \param[in] fp the file handle to read from.
 * \param[in] p a description of what is being read, used in the error message if the operation fails. */
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p) {
    if(fread(ptr,size,count,fp)!=count) {
        fprintf(stderr,"mesh: can't read %s from file\n",p);
        exit(1);
    }
}

/** Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void fatal_error(const char *p,int code) {
    fprintf(stderr,"Error: %s\n",p);
    exit(code);
}