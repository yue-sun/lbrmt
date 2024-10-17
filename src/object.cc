#include <cmath>
#include <cstdio>

#include "obj_field.hh"
#include "object.hh"

/** Calculate the initial velocity at a position, in terms of a constant
 * velocity plus an angular velocity around the object center.
 * \param[in] (X,Y) the reference map at the position.
 * \param[out] (ux,uy) the velocity. */
void obj_standard::velocity(double X,double Y,double &ux,double &uy) {
    ux=obj_ux-omega*Y;
    uy=obj_uy+omega*X;
}

/** Compute the level set function of the circle as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return the level set function. */
double obj_circle::phi(double X,double Y) {
    return sqrt(X*X+Y*Y)-cr;
}

/** Apply a gravitational acceleration to the circle.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phi the level set value at this position.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_circle::accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {
    accy-=grav;
}

/** Compute the level set function of the ellipse as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \paran[in] (sx_,sy_) the axis aspect ratio. 
 * \return the level set function. */
double obj_ellipse::phi(double X,double Y) {
    return sqrt(sx*X*X+sy*Y*Y)-cr;
}

/** Apply a gravitational acceleration to the ellipse.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phi the level set value at this position.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_ellipse::accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {
    accy-=grav;
}

/** Compute the level set function of the square as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return the level set function. */
double obj_square::phi(double X,double Y) {
    X=fabs(X)-cr;
    Y=fabs(Y)-cr;
    return X<0||Y<0?std::max(X,Y):sqrt(X*X+Y*Y);
}

/** Apply a gravitational acceleration to the square.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phi the level set value at this position.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_square::accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {
    accy-=grav;
}

/** Compute the level set function of the triangle as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return the level set function. */
double obj_triangle::phi(double X,double Y) {
    return sp(Y,-X,sp(X,Y,0)>0?(sp(X,Y,1)>0?0:2):(sp(X,Y,2)>0?1:0))-cr;
}

/** Apply a gravitational acceleration to the triangle.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phi the level set value at this position.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_triangle::accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {
    accy-=grav;
}

/** Compute the level set function of the flapper as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return the level set function. */
double obj_flapper::phi(double X,double Y) {
    double fX=fabs(X)-cl;
    return (fX<0?fabs(Y):sqrt(fX*fX+Y*Y))-cr;
}

/** Compute the current contraction amount of the flapper.
 * \param[in] time the current simulation time. */
void obj_flapper::pre_stress_setup(double time) {
    double a=sin(omega*time);
    a*=a*a;a*=a;a*=a;
    con=-amp*a;
}

/** Modify the deformation gradient tensor to apply the flapper swimming motion.
 * This default function applies no modification, but it can be overridden.
 * \param[in] (X,Y) the reference map at the position.
 * \param[in,out] F=[Xx Xy,Yx Yy] the deformation gradient tensor to modify. */
void obj_flapper::actuate(double X,double Y,double &Xx,double &Xy,double &Yx,double &Yy) {
    double fX=fabs(X)-al,phi=(fX<0.?fabs(Y):sqrt(fX*fX+Y*Y))-ar;
    if(phi<op->epsilon) {
        double sf=1.-op->trans_func(phi);
        double fac=exp(con*Y*sf);
        Xx*=fac;Xy*=1./fac;
        Yx*=fac;Yy*=1./fac;
    }
}

/** Compute the level set function of the smooth rotor as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return the level set function. */
double obj_smooth_rotor::phi(double X,double Y) {
    Y=fabs(Y);
    if(sp(X,Y,150)>0) {
        double yy=fabs(sp(X,Y,210));
        X=sp(X,Y,120);
        Y=yy;
    }
    return X>l?sqrt((X-l)*(X-l)+Y*Y)-r
              :(X>q?Y-r:s-sqrt((X-q)*(X-q)+(Y-r-s)*(Y-r-s)));
}

/** Calculate the rotation of the pivot prior to a stress computation.
 * \param[in] time the current simulation time. */
void obj_smooth_rotor::pre_stress_setup(double time) {
    double theta=th_max*(1.-cos(time*omega));
    cth=cos(theta);sth=sin(theta);
    // Compute the threshold for quickly determining if a point is inside the pivot region
    piv_thresh=piv_r+op->epsilon;
    piv_thresh*=piv_thresh;
}

/** Applies a rotating anchor force.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phi the level set value at this position.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_smooth_rotor::accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {
    double delx=x-cx,dely=y-cy,rsq=delx*delx+dely*dely;
    if(rsq<piv_thresh) {
        double rx=(X-cx)*cth-(Y-cy)*sth;
        double ry=(X-cx)*sth+(Y-cy)*cth;
        double K=0.01;
        accx-=K*(delx-rx);
        accy-=K*(dely-ry);
    }
}

/** Compute the level set function of the smooth rotor as a function of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return the level set function. */
double obj_seven_point::phi(double X,double Y) {
    const double fac=2*M_PI/7.;
    double theta=fabs(atan2(Y,X)),R=sqrt(X*X+Y*Y);
    theta=fabs(theta-fac*int(theta/fac+0.5));
    X=R*cos(theta)-cr;
    Y=R*sin(theta);
    return X*cos(0.25*fac)-Y*sin(0.25*fac)>0?sqrt(X*X+Y*Y):Y*cos(0.25*fac)+X*sin(0.25*fac);
}

/** Calculate the rotation of the pivot prior to a stress computation.
 * \param[in] time the current simulation time. */
void obj_seven_point::pre_stress_setup(double time) {
    double theta=th_max*(1-cos(time*omega));
    cth=cos(theta);sth=sin(theta);
    // Compute the threshold for quickly determining if a point is inside the pivot region
    piv_thresh=piv_r+op->epsilon;
    piv_thresh*=piv_thresh;
}

/** Applies a rotating anchor force.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phi the level set value at this position.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_seven_point::accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {
    double delx=x-cx,dely=y-cy,rsq=delx*delx+dely*dely;
    if(rsq<piv_thresh) {
        double rx=(X-cx)*cth-(Y-cy)*sth;
        double ry=(X-cx)*sth+(Y-cy)*cth;
        double K=0.01;
        accx-=K*(delx-rx);
        accy-=K*(dely-ry);
    }
}