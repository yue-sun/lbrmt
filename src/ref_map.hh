#ifndef REF_MAP
#define REF_MAP

#include <cmath>
#include <cstdio>

class ref_map {
    public:
        /** Reference map components. */
        double X,Y;
        /** Change of reference map components. */
        double cX,cY;
        /** Level set function value. */
        double phi;
        /** Extrapolation layer index, also used to label solid or fluid. */
        int cc;
        /** Blur zone half-width. */
        double epsilon;
        /** Object ID. */
        int id;
        /** Determinant of deformation gradient, det(F), i.e. solid volume. */
        double detF,invdetF;
        /** Extrapolation reset flag. */
        int extrap_reset;

        /** Instantiate class object. */
        ref_map() {}
        ref_map(double X_,double Y_,double phi_,int id_,double epsilon_) :
            X(X_), Y(Y_), cX(0.), cY(0.), phi(phi_), cc(phi<0.?0:-1),
            epsilon(epsilon_), id(id_), detF(1.), extrap_reset(0),
            tf1(0.5/epsilon), tf2(0.5/M_PI), tf3(M_PI/epsilon), lambda(trans_func(phi)) {}

        /** Destroy function. */
        virtual ~ref_map() {}

        /** Update reference map components. */
        inline void update() {
            X-=cX;
            Y-=cY;
        }

        /** Compute the solid fraction. */
        inline void solid_fraction() {
            lambda=trans_func(phi);
        }

        /** Calculate the smoothed Heaviside transition function.
         * 0: fully solid, 1: fully fluid.
         * \param[in] x the function argument.
         * \return the function value. */
        inline double trans_func(double x) {
            return x>epsilon?0:trans_func_in(x);
        }

        /** Calculate the smoothed Heaviside transition function,
         * assuming that the given argument corresponds to being inside the solid or the blur zone.
         * \param[in] x the function argument.
         * \return the function value. */
        inline double trans_func_in(double x) {
            return x<-epsilon?1:0.5-x*tf1-tf2*sin(tf3*x);
            
        }

    private:
        /** Constants that appear in the transition function calculations. */
        double tf1,tf2,tf3;

    public:
        /** Solid fraction. */
        double lambda;
};

#endif