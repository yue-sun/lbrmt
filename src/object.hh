#ifndef OBJECT_HH
#define OBJECT_HH

#include <cmath>
#include <cstdio>

class obj_field;
class lbrmt_2d;

/** Base class to define an solid within the simulation. */
class object {
    public:
        /** Destroy function. */
        virtual ~object() {}

        /** Level set function to define specific geometry.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual double phi(double X,double Y) = 0;

        /** Calculate the initial reference map at a position.
         * This default function uses the identity mapping, but it can be overridden
         * to apply different transforms.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            X=x;Y=y;
        }

        /** Calculate the reference map in the global space.
         * This default function uses the identity mapping, but it can be overridden
         * to apply different transforms.
         * \param[in] (X,Y) the coordinates of the reference map.
         * \param[out] (xx,yy) the coordinates of the reference map in the global space. */
        virtual void inverse_transform(double X,double Y,double &xx,double &yy) {
            xx=X;yy=Y;
        }


        /** Calculate the initial velocity at a position.
         * This default function sets the velocity to be zero, but it can be overridden.
         * \param[in] (X,Y) the reference map at the position.
         * \param[out] (ux,uy) the velocity. */
        virtual void velocity(double X,double Y,double &ux,double &uy) {
            ux=uy=0.;
        }

        /** Modify the deformation gradient tensor to apply actuation.
         * This default function applies no modification, but it can be overridden.
         * \param[in] (X,Y) the reference map at the position.
         * \param[in,out] F=[Xx Xy,Yx Yy] the deformation gradient tensor to modify. */
        virtual void actuate(double X,double Y,double &Xx,double &Xy,double &Yx,double &Yy) {}

        /** Link to the parent obj_field class to this object and perform additional initializations.
         * \param[in] lbrmt2d a reference to the parent lbrmt_2d class.
         * \param[in] op_ a pointer to the obj_field class representing this object. */
        inline void start(lbrmt_2d &lbrmt2d,obj_field *op_) {
            op=op_;
        }

        /** Perform any object-related output when a simulation frame is being stored.
         * This default function does nothing, but it can be overridden.
         * \param[in] k the frame number being stored.
         * \param[in] time the current simulation time. */
        virtual void write(int k,double time) {}

        /** Apply extra accelerations to the object, such as due to an anchor.
         * This default function does nothing, but it can be overriden.
         * \param[in] (x,y) the current position.
         * \param[in] (X,Y) the reference map at this position.
         * \param[in] phi the level set value at this position.
         * \param[in,out] (accx,accy) the acceleration vector to add to. */
        virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy) {}

        /** Compute the current contration amount for trig-related math functions.
         * This default function does nothing, but it can be overridden.
         * \param[in] time the current simulation time. */
        virtual void pre_stress_setup(double time) {}

    protected:
        /** A pointer to the corresponding obj_field class. */
        obj_field *op;
};

/** A special object type representing a solid filling the entire domain. */
struct obj_full : public object {
    /** Destroy function. */
    virtual ~obj_full() {}
    /** Level set function to define specific geometry.
     * \return (X,Y) the coordinates of the reference map. */
    virtual double phi(double X,double Y) {return -1;}
};

/** A base class for many objects, allowing them to be translated to different positions
 * and initialized with velocity and spin. */
struct obj_standard : public object {
    /** The position of the object center. */
    const double cx,cy;
    /** A length scale for the object size. */
    const double cr;
    /** The initial velocity of the object. */
    const double obj_ux,obj_uy;
    /** The initial angular velocity of the object. */
    const double omega;

    /** Initialize the class constants, setting the velocity parameters to zero.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale. */
    obj_standard(double cx_,double cy_,double cr_) :
        cx(cx_), cy(cy_), cr(cr_), obj_ux(0.), obj_uy(0.), omega(0.) {}

    /** Initialize the class constants, with nonzero velocity parameters.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale.
     * \param[in] (obj_ux_,obj_uy_) the initial object velocity.
     * \paran[in] omega_ the initial object angular velocity. */
    obj_standard(double cx_,double cy_,double cr_,double obj_ux_,double obj_uy_,double omega_) :
        cx(cx_), cy(cy_), cr(cr_), obj_ux(obj_ux_), obj_uy(obj_uy_), omega(omega_) {}

    /** Calculate the initial reference map at a position, applying a translation by the (cx,cy) vector.
     * \param[in] (x,y) the position to consider.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual void transform(double x,double y,double &X,double &Y) {
        X=x-cx;Y=y-cy;
    }

    /** Inverse the translation by the (cx,cy) vector in the initialization.
     * \param[in] (X,Y) the coordinates of the reference map.
     * \param[out] (xx,yy) the coordinates of the reference map in the global space. */
    virtual void inverse_transform(double X,double Y,double &xx,double &yy) {
        xx=X+cx;yy=Y+cy;
    }

    /** Calculate the initial velocity at a position.
     * \param[in] (X,Y) the reference map at the position.
     * \param[out] (ux,uy) the velocity. */
    virtual void velocity(double X,double Y,double &u,double &v);
};

/** A class to define a circle that experiences an acceleration in the y direction (e.g. due to gravity). */
struct obj_circle : public obj_standard {
    /** The acceleration to apply in the y direction. */
    const double grav;

    /** Initialize the class constants, setting the velocity parameters to zero.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale. 
     * \paran[in] grav_ the acceleration in the y direction. */
    obj_circle(double cx_,double cy_,double cr_,double grav_=0.) :
        obj_standard(cx_,cy_,cr_), grav(grav_) {}

    /** Initialize the class constants, with nonzero velocity parameters.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale.
     * \param[in] (obj_ux_,obj_uy_) the initial object velocity.
     * \paran[in] omega_ the initial object angular velocity.
     * \paran[in] grav_ the acceleration in the y direction. */
    obj_circle(double cx_,double cy_,double cr_,double obj_ux_,double obj_uy_,double omega_,double grav_=0.) :
        obj_standard(cx_,cy_,cr_,obj_ux_,obj_uy_,omega_), grav(grav_) {}

    /** Level set function to define specific geometry.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual double phi(double X,double Y);

    /** Apply extra accelerations to the object, such as due to an anchor.
     * \param[in] (x,y) the current position.
     * \param[in] (X,Y) the reference map at this position.
     * \param[in] phi the level set value at this position.
     * \param[in,out] (accx,accy) the acceleration vector to add to. */
    virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy);
};

/** A class to define an ellipse that experiences an acceleration in the y direction (e.g. due to gravity). */
struct obj_ellipse : public obj_standard {
    /** The aspect ratio of the object axis. */
    const double sx,sy;
    /** The acceleration to apply in the y direction. */
    const double grav;

    /** Initialize the class constants, setting the velocity parameters to zero.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale. 
     * \paran[in] (sx_,sy_) the axis aspect ratio. 
     * \paran[in] grav_ the acceleration in the y direction. */
    obj_ellipse(double cx_,double cy_,double cr_,double sx_,double sy_,double grav_=0.) :
        obj_standard(cx_,cy_,cr_), sx(sx_), sy(sy_), grav(grav_) {}

    /** Initialize the class constants, with nonzero velocity parameters.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale.
     * \paran[in] (sx_,sy_) the axis aspect ratio. 
     * \param[in] (obj_ux_,obj_uy_) the initial object velocity.
     * \paran[in] omega_ the initial object angular velocity.
     * \paran[in] grav_ the acceleration in the y direction. */
    obj_ellipse(double cx_,double cy_,double cr_,double sx_,double sy_,
        double obj_ux_,double obj_uy_,double omega_,double grav_=0.) :
        obj_standard(cx_,cy_,cr_,obj_ux_,obj_uy_,omega_), sx(sx_), sy(sy_), grav(grav_) {}

    /** Level set function to define specific geometry.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual double phi(double X,double Y);

    /** Apply extra accelerations to the object, such as due to an anchor.
     * \param[in] (x,y) the current position.
     * \param[in] (X,Y) the reference map at this position.
     * \param[in] phi the level set value at this position.
     * \param[in,out] (accx,accy) the acceleration vector to add to. */
    virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy);
};

/** A class to define a square that experiences an acceleration in the y direction (e.g. due to gravity). */
struct obj_square : public obj_standard {
    /** The acceleration to apply in the y direction. */
    const double grav;

    /** Initialize the class constants, setting the velocity parameters to zero.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale. 
     * \paran[in] grav_ the acceleration in the y direction. */
    obj_square(double cx_,double cy_,double cr_,double grav_=0.) :
        obj_standard(cx_,cy_,cr_), grav(grav_) {}

    /** Initialize the class constants, with nonzero velocity parameters.
     * \param[in] (cx_,cy_) the position of the object center.
     * \paran[in] cr_ the object length scale.
     * \param[in] (obj_ux_,obj_uy_) the initial object velocity.
     * \paran[in] omega_ the initial object angular velocity.
     * \paran[in] grav_ the acceleration in the y direction. */
    obj_square(double cx_,double cy_,double cr_,double obj_ux_,double obj_uy_,double omega_,double grav_=0.) :
        obj_standard(cx_,cy_,cr_,obj_ux_,obj_uy_,omega_), grav(grav_) {}

    /** Level set function to define specific geometry.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual double phi(double X,double Y);

    /** Apply extra accelerations to the object, such as due to an anchor.
     * \param[in] (x,y) the current position.
     * \param[in] (X,Y) the reference map at this position.
     * \param[in] phi the level set value at this position.
     * \param[in,out] (accx,accy) the acceleration vector to add to. */
    virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy);
};

/** A class to define a triangle that experiences an acceleration in the y direction (e.g. due to gravity). */
struct obj_triangle : public obj_standard {
    public:
        /** The acceleration to apply in the y direction. */
        const double grav;

        /** Initialize the class constants, setting the velocity parameters to zero.
         * \param[in] (cx_,cy_) the position of the object center.
         * \paran[in] cr_ the object length scale. 
         * \paran[in] theta the object angle. 
         * \paran[in] grav_ the acceleration in the y direction. */
        obj_triangle(double cx_,double cy_,double cr_,double theta,double grav_=0.) :
            obj_standard(cx_,cy_,cr_*0.5), grav(grav_) {
                for(int i=0;i<3;i++) {
                    nor[2*i]=cos(theta+2*M_PI/3*i);
                    nor[2*i+1]=sin(theta+2*M_PI/3*i);
                }
            }

        /** Level set function to define specific geometry.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual double phi(double X,double Y);

        /** Apply extra accelerations to the object, such as due to an anchor.
         * \param[in] (x,y) the current position.
         * \param[in] (X,Y) the reference map at this position.
         * \param[in] phi the level set value at this position.
         * \param[in,out] (accx,accy) the acceleration vector to add to. */
        virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy);

    private:
        /** THe normal vectors for the sides of triangle. */
        double nor[6];

        /** Computes the scalar product beteen the reference map and one of the normal vectors.
         * \param[in] (X,Y) the reference map.
         * \param[in] k the normal vector to use.
         * \return the scalr product. */
        inline double sp(double X,double Y,int k) {
            return X*nor[2*k]+Y*nor[2*k+1];
        }
};

/** A class to define a filament that can swim with a jellyfish-like motion. */
struct obj_flapper : public object {
    public:
        /** The x position of the flapper center. */
        const double cx;
        /** The y position of the flapper center. */
        const double cy;
        /** The radius of the flapper end caps. */
        const double cr;
        /** The half-length of the flapper (not including the end caps). */
        const double cl;
        /** The radius of the end caps of the actuated region. */
        const double ar;
        /** The half-length of the actuated region. */
        const double al;
        /** The maximum flapping amplitude. */
        const double amp;
        /** The angular frequency of the flapping motion. */
        const double omega;

        obj_flapper(double cx_,double cy_,double cr_,double cl_,double ar_,double al_,
                    double theta,double m_stretch,double T) :
            cx(cx_), cy(cy_), cr(cr_), cl(cl_), ar(ar_), al(al_), amp(log(m_stretch)/ar), omega(2.*M_PI/T),
            cth(cos(theta)), sth(sin(theta)) {}

        /** Calculate the initial reference map at a position, applying a translation by the (cx,cy) vector,
         * and applying a rotation.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            x-=cx;y-=cy;
            X=x*cth+y*sth;
            Y=-x*sth+y*cth;
        }

        /** Inverse the translation by the (cx,cy) vector and rotation in the initialization.
         * \param[in] (X,Y) the coordinates of the reference map.
         * \param[out] (xx,yy) the coordinates of the reference map in the global space. */
        virtual void inverse_transform(double X,double Y,double &xx,double &yy) {
            xx=X*cth-Y*sth+cx;yy=X*sth+Y*cth+cy;
        }

        /** Level set function to define specific geometry.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual double phi(double X,double Y);

        /** Compute the current contration amount for trig-related math functions.
         * This default function does nothing, but it can be overridden.
         * \param[in] time the current simulation time. */
        virtual void pre_stress_setup(double time);

        /** Modify the deformation gradient tensor to apply actuation.
         * This default function applies no modification, but it can be overridden.
         * \param[in] (X,Y) the reference map at the position.
         * \param[in,out] F=[Xx Xy,Yx Yy] the deformation gradient tensor to modify. */
        virtual void actuate(double X,double Y,double &Xx,double &Xy,double &Yx,double &Yy);
    private:
        /** Cosine of the flapper rotation angle. */
        const double cth;
        /** Sine of the flapper rotation angle. */
        const double sth;
        /** The current contraction level of the flapper. */
        double con;
};

/** A class to define a smooth rotor. */
struct obj_smooth_rotor : public object {
    public:
        /** The x position of the rotor center. */
        const double cx;
        /** The y position of the rotor center. */
        const double cy;
        /** The radius of the end caps. */
        const double r;
        /** The length of each prong from the origin. */
        const double l;
        /** The radius of the smoothed section between each prong. */
        const double s;
        /** The radius over which to apply the pivoting force. */
        const double piv_r;
        /** The rotation angle of the rotor. */
        const double rot;
        /** Half of the maximum theta value to spin the rotor to. */
        const double th_max;
        /** The angular velocity of the applied spinning term. */
        const double omega;

        obj_smooth_rotor(double cx_,double cy_,double r_,double l_,double s_,
            double piv_r_,double rot_,double th_max_,double omega_) : 
            cx(cx_), cy(cy_), r(r_), l(l_), s(s_), piv_r(piv_r_), rot(rot_),
            th_max(th_max_), omega(omega_), crot(cos(rot)), srot(sin(rot)), q((r+s)*sqrt(1./3)) {}

        /** Calculate the initial reference map at a position, applying a translation by the (cx,cy) vector,
         * and applying a rotation.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            x-=cx;y-=cy;
            X=x*crot+y*srot;
            Y=-x*srot+y*crot;
        }

        /** Inverse the translation by the (cx,cy) vector and rotation in the initialization.
         * \param[in] (X,Y) the coordinates of the reference map.
         * \param[out] (xx,yy) the coordinates of the reference map in the global space. */
        virtual void inverse_transform(double X,double Y,double &xx,double &yy) {
            xx=X*cth-Y*sth+cx;yy=X*sth+Y*cth+cy;
        }

        /** Level set function to define specific geometry.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual double phi(double X,double Y);

        /** Compute the current contration amount for trig-related math functions.
         * This default function does nothing, but it can be overridden.
         * \param[in] time the current simulation time. */
        virtual void pre_stress_setup(double time);

        /** Apply extra accelerations to the object, such as due to an anchor.
         * \param[in] (x,y) the current position.
         * \param[in] (X,Y) the reference map at this position.
         * \param[in] phi the level set value at this position.
         * \param[in,out] (accx,accy) the acceleration vector to add to. */
        virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy);
    private:
        /** A threshold used for rapid computations of the pivot region. */
        double piv_thresh;
        /** The cosine of the current angle of the rotor. */
        double cth;
        /** The sine of the current angle of the rotor. */
        double sth;
        /** The cosine of the rotation angle of the rotor. */
        double crot;
        /** The sine of the rotation angle of the rotor. */
        double srot;
        /** The cutoff distance for the smoothed section. */
        double q;
        /** Computes the scalar product between the reference map and one of
         * the normal vectors.
         * \param[in] (X,Y) the reference map.
         * \param[in] k the normal vector to use.
         * \return The scalar product. */
        inline double sp(double X,double Y,double theta) {
            double rad=M_PI/180.*theta;
            return X*cos(rad)+Y*sin(rad);
        }
};

/** A class to define a seven-pointed star. */
struct obj_seven_point : public object {
    public:
        /** The x position of the object center. */
        const double cx;
        /** The y position of the object center. */
        const double cy;
        /** The distance from the center to the star tips. */
        const double cr;
        /** The radius over which to apply the pivoting force. */
        const double piv_r;
        /** Half of the maximum theta value to spin the star to. */
        const double th_max;
        /** The angular velocity of the applied spinning term. */
        const double omega;

        obj_seven_point(double cx_,double cy_,double cr_,double piv_r_,
            double th_max_,double omega_) :
            cx(cx_), cy(cy_), cr(cr_), piv_r(piv_r_), th_max(th_max_), omega(omega_) {}

        /** Calculate the initial reference map at a position, applying a translation by the (cx,cy) vector,
         * and applying a rotation.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            x-=cx;y-=cy;
            X=x*cth+y*sth;
            Y=-x*sth+y*cth;
        }

        /** Inverse the translation by the (cx,cy) vector and rotation in the initialization.
         * \param[in] (X,Y) the coordinates of the reference map.
         * \param[out] (xx,yy) the coordinates of the reference map in the global space. */
        virtual void inverse_transform(double X,double Y,double &xx,double &yy) {
            xx=X*cth-Y*sth+cx;yy=X*sth+Y*cth+cy;
        }

        /** Level set function to define specific geometry.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual double phi(double X,double Y);

        /** Compute the current contration amount for trig-related math functions.
         * This default function does nothing, but it can be overridden.
         * \param[in] time the current simulation time. */
        virtual void pre_stress_setup(double time);

        /** Apply extra accelerations to the object, such as due to an anchor.
         * \param[in] (x,y) the current position.
         * \param[in] (X,Y) the reference map at this position.
         * \param[in] phi the level set value at this position.
         * \param[in,out] (accx,accy) the acceleration vector to add to. */
        virtual void accel(double x,double y,double X,double Y,double phi,double &accx,double &accy);
    private:
        /** A threshold used for rapid computations of the pivot region. */
        double piv_thresh;
        /** The cosine of the current angle of the spinner. */
        double cth;
        /** The sine of the current angle of the spinner. */
        double sth;
};

#endif