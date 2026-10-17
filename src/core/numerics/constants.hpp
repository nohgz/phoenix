/**
 * @file constants.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains all useful CGS constants and various pre-computed values.
 *        They're left as doubles and cast to the correct type at compile-time
 *        wherever used.
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 *
 */
#pragma once

#include <cmath>

namespace core {
  namespace constants {
    /**
     * @brief Proton mass \f$m_p\f$
     */
    inline constexpr double M_PROTON = 1.6726219e-24;
    /**
     * @brief Electron mass \f$m_e\f$
     */
    inline constexpr double M_ELECTRON = 9.1093829e-28;
    /**
     * @brief Electron charge \f$q_e\f$
     */
    inline constexpr double Q_ELECTRON = 4.80320425e-10;
    /**
     * @brief Boltzmann's constant, \f$k\f$
     */
    inline constexpr double K_BOLTZMANN = 1.3806488e-16;
    /**
     * @brief The speed of light, \f$c\f$
     */
    inline constexpr double SPEED_OF_LIGHT = 2.99792458e10;
    /**
     * @brief Planck's constant, \f$h\f$
     */
    inline constexpr double H_PLANCK = 6.62606885e-27;
    /**
     * @brief Gravitational constant, \f$G\f$
     */
    inline constexpr double G_GRAV = 6.674e-8;
    /**
     * @brief Solar mass in g
     */
    inline constexpr double M_SUN = 1.989e33;
    /**
     * @brief Kiloparsec in cm
     */
    inline constexpr double KPC = 3.086e21;
    /**
     * @brief erg/s/cm^2/Hz -> Jy
     */
    inline constexpr double JANSKY = 1e23;


    // ******* PRE COMPUTED CONSTANTS *******
    /**
     * @brief \f$h / (m_e c^2)\f$: dimensionless photon energy per Hz
     */
    inline constexpr double H_OVER_MEC2 = H_PLANCK / (M_ELECTRON * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    /**
     * @brief \f$3 q_e / (4\pi m_e c)\f$
     */
    inline constexpr double NU_CRIT_COEFF = 3. * Q_ELECTRON / (4. * M_PI * M_ELECTRON * SPEED_OF_LIGHT);
    /**
     * @brief \f$q_e^2 / (2\sqrt{3} c)\f$: emissivity prefactor
     */
    inline constexpr double JNU_COEFF = Q_ELECTRON * Q_ELECTRON / (2. * 1.73205080756887729353 * SPEED_OF_LIGHT);
    /**
     * @brief \f$4\pi q_e^2 / m_e\f$: plasma frequency squared per \f$n_e\f$
     */
    inline constexpr double PLASMA_FREQ2_COEFF = 4. * M_PI * Q_ELECTRON * Q_ELECTRON / M_ELECTRON;
    /**
     * @brief \f$q_e / (m_e c)\f$: electron gyro frequency per \f$B\f$
     */
    inline constexpr double GYRO_FREQ_COEFF = Q_ELECTRON / (M_ELECTRON * SPEED_OF_LIGHT);
  }
}