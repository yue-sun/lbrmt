/** This C++ code is adapted from Chris Rycroft VDropImpact sample code to
 * parse simulation config file into the solver.
 * https://github.com/chr1shr/vdropimpact/blob/master/fileinfo.cc
 * https://github.com/chr1shr/vdropimpact/blob/master/fileinfo.hh */

#include "convert_cfg.hh"
#include "file.hh"

/** Load .cfg file and convert input texts into simulation parameters.
 * \param[in] infile input .cfg file.
 * \param[in,out] outdir output directory name.
 * \param[in,out] outdir_mmap output multimaps directory name.
 * \param[in,out] niters number of total simulation iterations.
 * \param[in,out] nout number of output frames.
 * \param[in,out] nx number of grid points in x direction.
 * \param[in,out] ny number of grid points in y direction.
 * \param[in,out] epsilon blur zone half-width.
 * \param[in,out] bctype fluid boundary conditions.
 * \param[in,out] rho_f dimensionless fluid density.
 * \param[in,out] Re Reynolds number.
 * \param[in,out] tau LB relaxation time.
 * \param[in,out] nu_ph physical fluid kinematic viscosity.
 * \param[in,out] ux0 initial fluid velocity in x direction.
 * \param[in,out] uy0 initial fluid velocity in y direction.
 * \param[in,out] fx_fluid_ph physical fluid force density in x direction.
 * \param[in,out] fy_fluid_ph physical fluid force density in y direction.
 * \param[in,out] dx physical grid spacing.
 * \param[in,out] dt physical time step.
 * \param[in,out] C_rho density scale.
 * \param[in,out] nobjs number of solid objects. */
void convert_cfg(const char* infile,char* &outdir,char* &outdir_mmap,
    int &niters,int &nout,int &nx,int &ny,double &epsilon,int &bctype,
    double &rho_f,double &Re,double &tau,double &nu_ph,
    double &ux0,double &uy0,double &fx_fluid_ph,double &fy_fluid_ph,
    double &dx,double &dt,double &C_rho,int &nobjs) {
    // Check that the input filename ends in '.cfg'
    int l=strlen(infile),ln=1;
    if(l<4) fatal_error("Filename is too short",1);
    const char* ip=infile+l-4;
    if(*ip!='.'||ip[1]!='c'||ip[2]!='f'||ip[3]!='g')
        fatal_error("Filename must end in '.cfg'",1);

    // Assemble output filename by replacing '.cfg' with '.out'
    outdir=new char[l+1];
    memcpy(outdir,infile,l-3);
    char *fp=outdir+l-3;
    *fp='o';fp[1]='u';fp[2]='t';fp[3]=0;

    // Assemble output multimaps filename by appending '/mmap'
    int nl=strlen(outdir);
    outdir_mmap=new char[nl+6];
    memcpy(outdir_mmap,outdir,nl);
    char *fp_mmap=outdir_mmap+l;
    *fp_mmap='/';fp_mmap[1]='m';fp_mmap[2]='m';
    fp_mmap[3]='a';fp_mmap[4]='p';fp_mmap[5]=0;
    
    // Open the input file and read
    FILE *f=safe_fopen(infile,"r");
    char *buf=new char[cfg_buf_size],*bp;
    while(!feof(f)) {
        if(fgets(buf,cfg_buf_size,f)==NULL) break;

        // Locate comments and remove by replacing comment character
        // with a null character
        bp=buf;
        while((*bp)!=0) {
            if(*bp=='#') {*bp=0;break;}
            bp++;
        }

        // Separate entries in the file by tabs and spaces; if no entries are
        // available then skip this line
        bp=strtok(buf," \t\n");
        if(bp==NULL) {ln++;continue;}

        // Look for a known keyword, and read in any extra values
        if(se(bp,"niters")) niters=final_int(ln);
        else if(se(bp,"nout")) nout=final_int(ln);
        else if(se(bp,"nx")) {
            nx=final_int(ln);
            if(nx<=0) fatal_error("Grid dimensions nx must be positive",1);
        } else if(se(bp,"ny")) {
            ny=final_int(ln);
            if(ny<=0) fatal_error("Grid dimensions ny must be positive",1);
        } else if(se(bp,"epsilon")) epsilon=final_double(ln);
        else if(se(bp,"bctype")) bctype=final_int(ln);
        else if(se(bp,"rho_f")) rho_f=final_double(ln);
        else if(se(bp,"Re")) Re=final_double(ln);
        else if(se(bp,"tau")) tau=final_double(ln);
        else if(se(bp,"nu_ph")) nu_ph=final_double(ln);
        else if(se(bp,"ux0")) ux0=final_double(ln);
        else if(se(bp,"uy0")) uy0=final_double(ln);
        else if(se(bp,"fx_fluid_ph")) fx_fluid_ph=final_double(ln);
        else if(se(bp,"fy_fluid_ph")) fy_fluid_ph=final_double(ln);
        else if(se(bp,"dx")) dx=final_double(ln);
        else if(se(bp,"dt")) dt=final_double(ln);
        else if(se(bp,"C_rho")) C_rho=final_double(ln);
        else if(se(bp,"nobjs")) nobjs=final_int(ln);
        ln++;
    }
    delete [] buf;
    fclose(f);

    // Check that all parameters have been specified
    check_invalid(niters,"niters");
    check_invalid(nout,"nout");
    if(nx==-1) fatal_error("Grid dimensions nx not set",1);
    if(ny==-1) fatal_error("Grid dimensions ny not set",1);
    check_invalid(epsilon,"epsilon");
    check_invalid(bctype,"bctype");
    check_invalid(rho_f,"rho_f");
    check_invalid(Re,"Re");
    check_invalid(tau,"tau");
    check_invalid(nu_ph,"nu_ph");
    check_invalid(ux0,"ux0");
    check_invalid(uy0,"uy0");
    check_invalid(fx_fluid_ph,"fx_fluid_ph");
    check_invalid(fy_fluid_ph,"fy_fluid_ph");
    check_invalid(dx,"dx");
    check_invalid(dt,"dt");
    check_invalid(C_rho,"C_rho");
    check_invalid(nobjs,"nobjs");
}

/** Load .cfg file and convert input texts into simulation parameters
 * for solid objects.
 * \param[in] infile input .cfg file.
 * \param[in,out] obj_type solid object type.
 * \param[in,out] rho_s dimensionless solid density.
 * \param[in,out] G dimensionless solid shear modulus.
 * \param[in,out] cx dimensionless center of solid object in x direction.
 * \param[in,out] cy dimensionless center of solid object in y direction.
 * \param[in,out] cr dimensionless size of solid object. */
void convert_solid_cfg(const char* infile,int &obj_type,double &rho_s,
    double &G,double &cx,double &cy,double &cr) {
    // Check that the input filename ends in '.cfg'
    int l=strlen(infile),ln=1;
    if(l<4) fatal_error("Filename is too short",1);
    const char* ip=infile+l-4;
    if(*ip!='.'||ip[1]!='c'||ip[2]!='f'||ip[3]!='g')
        fatal_error("Filename must end in '.cfg'",1);
    
    // Open the input file and read
    FILE *f=safe_fopen(infile,"r");
    char *buf=new char[cfg_buf_size],*bp;
    while(!feof(f)) {
        if(fgets(buf,cfg_buf_size,f)==NULL) break;

        // Locate comments and remove by replacing comment character
        // with a null character
        bp=buf;
        while((*bp)!=0) {
            if(*bp=='#') {*bp=0;break;}
            bp++;
        }

        // Separate entries in the file by tabs and spaces; if no entries are
        // available then skip this line
        bp=strtok(buf," \t\n");
        if(bp==NULL) {ln++;continue;}

        // Look for a known keyword, and read in any extra values
        if(se(bp,"obj_type")) obj_type=final_int(ln);
        else if(se(bp,"rho_s")) rho_s=final_double(ln);
        else if(se(bp,"G")) G=final_double(ln);
        else if(se(bp,"cx")) cx=final_double(ln);
        else if(se(bp,"cy")) cy=final_double(ln);
        else if(se(bp,"cr")) cr=final_double(ln);
        ln++;
    }
    delete [] buf;
    fclose(f);

    // Check that all parameters have been specified
    check_invalid(obj_type,"obj_type");
    check_invalid(rho_s,"rho_s");
    check_invalid(G,"G");
    check_invalid(cx,"cx");
    check_invalid(cy,"cy");
    check_invalid(cr,"cr");
}

/** Find the next token in a string and interpret it as a double precision 
 * floating point number. If none is availble, it gives an error message.
 * \param[in] ln the current line number. */
inline double next_double(int ln) {
    return atof(next_token(ln));
}

/** Find the next token in a string, interpret it as a double precision 
 * floating point number, and checks that there are no subsequent values.
 * \param[in] ln the current line number. */
inline double final_double(int ln) {
    double temp=next_double(ln);
    check_no_more(ln);
    return temp;
}

/** Find the next token in a string, interpret it as an integer, and check
 * that there are no subsequent values.
 * \param[in] ln the current line number. */
inline int final_int(int ln) {
    int temp=atoi(next_token(ln));
    check_no_more(ln);
    return temp;
}

/** Test to see if two strings are equal.
 * \param[in] p1 a pointer to the first string.
 * \param[in] p2 a pointer to the second string.
 * \return True if they are equal, false otherwise. */
inline bool se(const char *p1,const char *p2) {
    return strcmp(p1,p2)==0;
}

/** Find the next token in a string, and if none is availble, raise error.
 * \param[in] ln the current line number. */
char* next_token(int ln) {
    char *temp=strtok(NULL," \t\n");
    if(temp==NULL) {
        fprintf(stderr,"Not enough arguments at input line %d\n",ln);
        exit(1);
    }
    return temp;
}

/** Check that there are no subsequent values.
 * \param[in] ln the current line number. */
void check_no_more(int ln) {
    if(strtok(NULL," \t\n")!=NULL) {
        fprintf(stderr,"Too many arguments at input line %d\n",ln);
        exit(1);
    }
}

/** Check that a parameter is a valid positive value.
 * \param[in] val the value to check.
 * \param[in] p the name of the value. */
void check_invalid(double val,const char *p) {
    if(val<0) {
        fprintf(stderr,"Value of %s either invalid or not set\n",p);
        exit(1);
    }
}

/** Generate a random double between two doubles.
 * \param[in] a,b two doubles.
 * \param[out] rand the random double between a and b. */
double rand_ab(double a,double b) {
    double random=((double) rand())/(double) RAND_MAX;
    double diff=b-a;
    double rand=a+random*diff;
    return rand;
}