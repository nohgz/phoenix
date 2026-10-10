/**
 * @file metric.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief The MKSBHAC (and todo others) metrics with tensor helpers
 * @version 1.0
 * @date 2026-10-09
 *
 * GNU GPL v3.0
 *
 */

#pragma once

#include "core/numerics/headers.hpp"

namespace core {
  namespace metric {
    /**
     * @brief Covariant MKSBHAC metric at X_u = (t, ln r, x2, phi).
     *
     * @param a      Black hole spin
     * @param hslope MKSBHAC theta-stretching parameter
     * @param X_u    Position
     * @param g_dd   Output covariant metric
     */
    template <class T>
    GRRT_SHARED_FUNCTION void metric_dd(
      T a,
      T hs,
      const T X_u[4],
      T g_dd[4][4]
    ) {
      const T two = T(2), one = T(1);
      const T r = std::exp(X_u[1]);
      T s2x, c2x, s, c;
      sincos_t(two * X_u[2], &s2x, &c2x);
      const T th = X_u[2] + T(0.5) * hs * s2x;
      sincos_t(th, &s, &c);
      const T s2 = s * s;
      const T rho2 = r * r + a * a * c * c;
      const T hfac = one + hs * c2x;
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          g_dd[i][j] = T(0);
        }
      }
      g_dd[0][0] = -one + two * r / rho2;
      g_dd[0][1] = g_dd[1][0] = (two * r / rho2) * r;
      g_dd[0][3] = g_dd[3][0] = -two * a * r * s2 / rho2;
      g_dd[1][1] = (one + two * r / rho2) * r * r;
      g_dd[1][3] = g_dd[3][1] = -a * s2 * (one + two * r / rho2) * r;
      g_dd[2][2] = rho2 * hfac * hfac;
      g_dd[3][3] = s2 * (rho2 + a * a * s2 * (one + two * r / rho2));
    }

    /**
     * @brief Contravariant MKSBHAC metric at X_u = (t, ln r, x2, phi).
     *
     * @param a    Black hole spin
     * @param hs   MKSBHAC theta-stretching parameter (aka hslope)
     * @param X_u  Position
     * @param g_uu Output contravariant metric
     */
    template <class T>
    GRRT_SHARED_FUNCTION void metric_uu(
      T a,
      T hs,
      const T X_u[4],
      T g_uu[4][4]
    ) {
      const T two = T(2), one = T(1);
      const T r = std::exp(X_u[1]);
      T s2x, c2x, s, c;
      sincos_t(two * X_u[2], &s2x, &c2x);
      const T th = X_u[2] + T(0.5) * hs * s2x;
      sincos_t(th, &s, &c);
      const T s2 = s * s;
      const T rho2 = r * r + a * a * c * c;
      const T hfac = one + hs * c2x;
      const T irho2 = one / rho2;
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          g_uu[i][j] = T(0);
        }
      }
      g_uu[0][0] = -one - two * r * irho2;
      g_uu[0][1] = g_uu[1][0] = two * irho2;
      g_uu[1][1] = irho2 * (r * (r - two) + a * a) / (r * r);
      g_uu[1][3] = g_uu[3][1] = a * irho2 / r;
      g_uu[2][2] = irho2 / (hfac * hfac);
      g_uu[3][3] = irho2 / s2;
    }

    /**
     * @brief Covariant & Contravariant MKSBHAC metrics at X_u = (t, ln r, x2, phi).
     *
     * @param a    Black hole spin
     * @param hs   MKSBHAC theta-stretching parameter (aka hslope)
     * @param X_u  Position
     * @param g_dd Output covariant metric
     * @param g_uu Output contravariant metric
     */
    template <class T>
    GRRT_SHARED_FUNCTION void metric_dd_uu(
      T a,
      T hs,
      const T X_u[4],
      T g_dd[4][4],
      T g_uu[4][4]
    ) {
      const T two = T(2), one = T(1);
      const T r = std::exp(X_u[1]);
      T s2x, c2x, s, c;
      sincos_t(two * X_u[2], &s2x, &c2x);
      const T th = X_u[2] + T(0.5) * hs * s2x;
      sincos_t(th, &s, &c);
      const T s2 = s * s;
      const T rho2 = r * r + a * a * c * c;
      const T irho2 = one / rho2;
      const T hfac = one + hs * c2x;

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          g_dd[i][j] = T(0);
          g_uu[i][j] = T(0);
        }
      }

      g_dd[0][0] = -one + two * r / rho2;
      g_dd[0][1] = g_dd[1][0] = (two * r / rho2) * r;
      g_dd[0][3] = g_dd[3][0] = -two * a * r * s2 / rho2;
      g_dd[1][1] = (one + two * r / rho2) * r * r;
      g_dd[1][3] = g_dd[3][1] = -a * s2 * (one + two * r / rho2) * r;
      g_dd[2][2] = rho2 * hfac * hfac;
      g_dd[3][3] = s2 * (rho2 + a * a * s2 * (one + two * r / rho2));

      g_uu[0][0] = - one - two * r / rho2;
      g_uu[0][1] = g_uu[1][0] = two * irho2;
      g_uu[1][1] = irho2 * (r * (r - two) + a * a) / (r * r);
      g_uu[1][3] = g_uu[3][1] = a * irho2 / r;
      g_uu[2][2] = irho2 / (hfac * hfac);
      g_uu[3][3] = irho2 / s2;
    }
  }
}