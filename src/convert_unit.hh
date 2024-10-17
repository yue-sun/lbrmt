#ifndef CONVERT_UNIT_HH
#define CONVERT_UNIT_HH

/** Convert the shear modulus G from physical units to LB units.
 * \param[in] G_ph the physical shear modulus.
 * \param[in] C_G the conversion factor of shear modulus.
 * \return the LB shear modulus. */
inline double convert_G(double G_ph,double C_G) {
    return G_ph/C_G;
}

/** Convert general force density from physical units to LB units.
 * \param[in] force_ph the physical force density.
 * \param[in] C_g the conversion factor of force density.
 * \return the LB force density. */
inline double convert_force(double force_ph,double C_g) {
    return force_ph/C_g;
}

/** Convert sedimentation force density from physical units to LB units.
 * \param[in] rho_s the LB fluid density.
 * \param[in] rho_s the LB solid density.
 * \param[in] C_g the conversion factor of force density.
 * \return the LB sedimentation force density. */
inline double convert_set_force(double rho_f,double rho_s,double C_g) {
    // Sedimentation force is implemented by factoring out background fluid pressure
    double force_ph=9.81*(1.-rho_f/rho_s);
    return force_ph/C_g;
}

#endif