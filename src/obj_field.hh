#ifndef OBJ_FIELD_HH
#define OBJ_FIELD_HH

#include <cmath>
#include <vector>

#include "object.hh"

/** A structure containing material constants for a solid. */
struct mat_const {
    /** The shear modulus. */
    double G;
    /** The solid density. */
    double rho_s;
    /** Whether the solid is initialized to its own velocity, or inherits it from the fluid. */
    bool set_velocity;

    mat_const(double G_,double rho_s_,bool set_velocity_=true) :
        G(G_), rho_s(rho_s_), set_velocity(set_velocity_) {}
    /** Destroy function. */
    ~mat_const() {}
};

/** A class describing all data needed to simulate a single solid. */
class obj_field {
    public:
        /** The shear modulus. */
        const double G;
        /** The solid density. */
        const double rho_s;
        /** The blur zone half-width. */
        const double epsilon;
        /** Whether the solid is initialized to its own velocity, or inherits it from the fluid. */
        const bool set_velocity;

        /** A pointer to the object class containing details about the object's geometry and other attributes. */
        object* obj;
        /** Vectors for constructing layers in the extrapolation zone (swap between l and l2 for each layer). */
        std::vector<int> l,l2;

        /** Instantitate class object. */
        obj_field(object* obj_,mat_const *mc,double epsilon_) :
            G(mc->G), rho_s(mc->rho_s), epsilon(epsilon_), set_velocity(mc->set_velocity), obj(obj_),
            tf1(0.5/epsilon), tf2(0.5/M_PI), tf3(M_PI/epsilon) {}
        /** Destroy function. */
        ~obj_field() {}

        /** Calculate the smoothed Heaviside transition function.
         * 0: fully solid, 1: fully fluid.
         * \param[in] x the function argument.
         * \return the function value. */
        inline double trans_func(double x) {
            return x>epsilon?1.:trans_func_in(x);
        }

        /** Calculate the smoothed Heaviside transition function,
         * assuming that the given argument corresponds to being inside the solid or the blur zone.
         * \param[in] x the function argument.
         * \return the function value. */
        inline double trans_func_in(double x) {
            return x<-epsilon?0.:0.5+tf1*x+tf2*sin(tf3*x);
        }

    private:
        /** Constants that appear in the transition function calculations. */
        const double tf1,tf2,tf3;
};

#endif