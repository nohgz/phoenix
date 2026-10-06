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
  }
}