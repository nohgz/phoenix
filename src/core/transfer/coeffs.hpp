/**
 * @file coeffs.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Thermal synchrotron emission / absorption / Faraday coefficients.
 * @version 1.2
 * @date 2026-10-13
 *
 * GNU GPL v3.0
 */

#pragma once

#include "core/numerics/headers.hpp"

namespace core {
  namespace coeffs {
    /**
     * @brief Polarized transfer coefficients.
     *        Tetrad formalism lets us align them with the projected B,
     *        so jU = aU = rU = 0.
     */
    struct PolCoeffs {
      scalar jI, jQ, jU, jV;  // emission
      scalar aI, aQ, aU, aV;  // absorption
      scalar rQ, rU, rV;      // Faraday conversion (Q) / rotation (V)
    };

    /**
     * @brief Zero out the transfer coeffs for initialization
     */
    GRRT_SHARED_FORCEINLINE PolCoeffs polcoeffs_zero(void) {
      PolCoeffs c;
      c.jI = c.jQ = c.jU = c.jV = scalar(0);
      c.aI = c.aQ = c.aU = c.aV = scalar(0);
      c.rQ = c.rU = c.rV = scalar(0);
      return c;
    }

    // Compile time constants to avoid excess expensive runtime calcs
    // 2^(-1/3) / 3^(23/6) * 1e4 * pi
    constexpr double COEF_C0 = 369.69342869090815;
    constexpr double LOG120 = 4.787491742782046;

    GRRT_SHARED_FUNCTION scalar I_I(scalar x) {
      const scalar cbrtx = std::cbrt(x);
      const scalar icbrtx = scalar(1) / cbrtx;  // x^(1/3), x^(-1/3)
      return scalar(2.5651) * (scalar(1) + scalar(1.92) * icbrtx + scalar(0.9977) * icbrtx * icbrtx) * std::exp(scalar(-1.8899) * cbrtx);
    }

    GRRT_SHARED_FUNCTION scalar I_Q(scalar x) {
      const scalar cbrtx = std::cbrt(x);
      const scalar icbrtx = scalar(1) / cbrtx;
      return scalar(2.5651) * (scalar(1) + scalar(0.93193) * icbrtx + scalar(0.499873) * icbrtx * icbrtx) * std::exp(scalar(-1.8899) * cbrtx);
    }

    GRRT_SHARED_FUNCTION scalar I_V(scalar x) {
      const scalar cbrtx = std::cbrt(x);
      const scalar icbrtx = scalar(1) / cbrtx;
      return (scalar(1.81348) * icbrtx * icbrtx * icbrtx + scalar(3.42319) * icbrtx * icbrtx + scalar(0.0292545) / std::sqrt(x) +
              scalar(2.03773) * icbrtx) * std::exp(scalar(-1.8899) * cbrtx);
    }

    GRRT_SHARED_FUNCTION scalar f_m(scalar X) {
      const scalar lx = std::log(X);
      const scalar ex = std::exp(-X / scalar(47.2));
      return scalar(2.011) * std::exp(-std::exp(scalar(1.035) * lx) / scalar(4.7)) -
             std::cos(X * scalar(0.5)) * std::exp(-std::exp(scalar(1.2) * lx) / scalar(2.73)) - scalar(0.011) * ex +
             (scalar(0.011) * ex - scalar(COEF_C0) * std::exp(scalar(-8. / 3.) * lx)) * scalar(0.5) *
                 (scalar(1) + std::tanh(scalar(10) * (lx - scalar(LOG120))));
    }

    /**
     * @brief Planck function B_nu(nu, Theta_e).
     */
    GRRT_SHARED_FUNCTION scalar planck_nu(scalar nu, scalar te) {
      const scalar pre = scalar(2. * constants::H_PLANCK) * nu * nu * nu /
                               scalar(constants::SPEED_OF_LIGHT * constants::SPEED_OF_LIGHT);
      const scalar xx = scalar(constants::H_OVER_MEC2) * nu / te;
      return IS_FP32 ? pre / std::expm1(xx) : pre / (std::exp(xx) - scalar(1));
    }

    /**
     * @brief Raw cgs coefficients (j_I, j_Q, j_V, a_*, rho_Q, rho_V).
     *
     * Float note: (2 pi nu)^4 ~ 4e48 overflows float. The Faraday terms are
     * formed from the dimensionless ratio omega0 / (2 pi nu) instead of the
     * raw powers (identical in exact arithmetic).
     *
     * @param te        Electron temperature [m_e c^2]
     * @param ne        Electron density [cm^-3]
     * @param nu        Plasma-frame frequency [Hz]
     * @param B         |B| [G]
     * @param thB       Angle(k, B) in the plasma frame
     */
    GRRT_SHARED_FORCEINLINE PolCoeffs thermal_raw(
      scalar te,
      scalar ne,
      scalar nu,
      scalar B,
      scalar thB
    ) {
      PolCoeffs r = polcoeffs_zero();
      scalar sth, cth;
      sincos_t(thB, &sth, &cth);
      const scalar ite = scalar(1) / te;
      const scalar nu_c = scalar(constants::NU_CRIT_COEFF) * B * sth * te * te;
      const scalar x = nu / nu_c;
      const scalar pref = scalar(constants::JNU_COEFF) * ne * nu * ite * ite;  // emissivity prefactor / theta_e^2
      r.jI = pref * I_I(x);
      r.jQ = pref * I_Q(x);
      r.jV = scalar(4. / 3.) * pref * ite * (cth / sth) * I_V(x);  // 2 ne e^2 nu / (3 sqrt3 c) / (tan th theta_e^3)
      const scalar Bnu = planck_nu(nu, te);
      const scalar iB = scalar(1) / Bnu;
      r.aI = r.jI * iB;
      r.aQ = r.jQ * iB;
      r.aV = r.jV * iB;
      const scalar iw = scalar(1) / (scalar(2 * M_PI) * nu);  // 1 / w
      const scalar wp2 = scalar(constants::PLASMA_FREQ2_COEFF) * ne;
      const scalar omega0 = scalar(constants::GYRO_FREQ_COEFF) * B;
      const scalar Xe =
          te * std::sqrt(scalar(1.4142135623730951) * sth * (scalar(1.e3 / (2. * M_PI)) * omega0 / nu));
      // only the ratios K1/K2 and K0/K2 enter the Faraday coefficients
      scalar kr1, kr0;
      if (ite < scalar(1. / 5.)) {
        const scalar k0 = -std::log(ite / scalar(2)) - scalar(0.5772), k1 = scalar(1) / ite,
                           k2 = scalar(2) / ite / ite;
        kr1 = k1 / k2;
        kr0 = (k2 > scalar(0)) ? k0 / k2 : scalar(1);
      } else {
        bessel::bessel_ratios(ite, &kr1, &kr0);
      }
      const scalar q = omega0 * iw;  // omega0 / (2 pi nu)
      const scalar wpw = wp2 * iw;   // wp^2 / w
      r.rQ = (scalar(0.5) / scalar(constants::SPEED_OF_LIGHT)) * wpw * q * q * f_m(Xe) * (kr1 + scalar(6) * te) * sth * sth;
      const scalar shgm = scalar(1) - scalar(0.11) * std::log(scalar(1) + scalar(0.035) * Xe);
      r.rV = (scalar(1) / scalar(constants::SPEED_OF_LIGHT)) * wpw * q * (kr0 * shgm) * cth;
      return r;
    }

    /**
     * @brief Plasma-frame cgs coefficients -> Lorentz-invariant form
     *
     * @param c         Raw coefficients
     * @param nu_p      Plasma-frame frequency
     * @param jscale    Multiplies the emissivities
     */
    GRRT_SHARED_FORCEINLINE PolCoeffs invariant(PolCoeffs c, scalar nu_p, scalar jscale) {
      const scalar js = jscale / nu_p / nu_p;
      c.jI *= js;
      c.jQ *= js;
      c.jV *= js;
      c.aI *= nu_p;
      c.aQ *= nu_p;
      c.aV *= nu_p;
      c.rQ *= nu_p;
      c.rV *= nu_p;
      const scalar pf = std::sqrt(c.jQ * c.jQ + c.jV * c.jV) / c.jI;
      if (pf > scalar(1)) {
        c.jQ /= (pf + scalar(0.005));
        c.jU /= (pf + scalar(0.005));
        c.jV /= (pf + scalar(0.005));
      }
      return c;
    }

    /**
     * @brief Rotates a plasma-frame coefficient set into the camera polarization basis.
     *
     * @param c   Coefficients in the plasma basis
     * @param c2  cos 2psi (chi_from_F1)
     * @param s2  sin 2psi
     */
    GRRT_SHARED_FORCEINLINE PolCoeffs rotate_to_camera(const PolCoeffs* c, scalar c2, scalar s2) {
      PolCoeffs r = *c;
      r.jQ = c2 * c->jQ + s2 * c->jU;
      r.jU = -s2 * c->jQ + c2 * c->jU;
      r.aQ = c2 * c->aQ + s2 * c->aU;
      r.aU = -s2 * c->aQ + c2 * c->aU;
      r.rQ = c2 * c->rQ + s2 * c->rU;
      r.rU = -s2 * c->rQ + c2 * c->rU;
      return r;
    }
  }
}