#pragma once

#include "core/numerics/headers.hpp"

namespace core {
  namespace geodesic {
    /**
     * @brief hi + lo = (x[4], k[4], F1[4]) in FP32. lo goes away in FP64 as its
     *        for TwoSum in FP32 only to compensate for error
     */
    struct Geodesic {
      scalar hi[12];
      scalar lo[12];
    };

    /**
     * @brief One lightpath sample to pass from geodesic to the plasma stage.
     */
    struct Sample {
      int i;
      // Position vector
      scalar x[4];
      // Wave vector
      scalar k[4];
      // Polarization vector
      scalar F1[4];
      // Stepsize
      scalar dl;
    };

    /**
     * @brief Geometric quantities shared by all RHS evaluations at one point.
     */
    struct Geom {
      scalar a, a2;
      scalar r, r2, tp, tp2, tpp;
      scalar s, c, s2, s2th, cot;
      scalar Sig, Sb, Delta, A, B, C;
    };

    /**
     * @brief Initial state from the double camera ray start.
     *
     * @param st       Output state
     * @param photon_u x_u, k_u
     * @param F1_u     Camera polarization vector
     */
    GRRT_SHARED_FUNCTION void init(
      Geodesic* geo,
      const double photon_u[8],
      const double F1_u[4]
    ) {
      for (int i = 0; i < 8; i++) geo->hi[i] = scalar(photon_u[i]);
      for (int i = 0; i < 4; i++) geo->hi[8 + i] = scalar(F1_u[i]);
      for (int i = 0; i < 12; i++) geo->lo[i] = scalar(0);
      if (IS_FP32) {
        for (int i = 0; i < 8; i++) geo->lo[i] = scalar(photon_u[i] - (double)geo->hi[i]);
        for (int i = 0; i < 4; i++) geo->lo[8 + i] = scalar(F1_u[i] - (double)geo->hi[8 + i]);
      }
    }

    /**
     * @brief Fill the sample with the current geodesic information
     */
    GRRT_SHARED_FORCEINLINE void view(const Geodesic* geo, Sample* smp) {
      for (int i = 0; i < 4; i++) {
        smp->x[i]  = scalar(geo->hi[i]);
        smp->k[i]  = scalar(geo->hi[4 + i]);
        smp->F1[i] = scalar(geo->hi[8 + i]);
        if (IS_FP32) {
          smp->x[i]  += scalar(geo->lo[i]);
          smp->k[i]  += scalar(geo->lo[4 + i]);
          smp->F1[i] += scalar(geo->lo[8 + i]);
        }
      }
    }

    /**
     * @brief Populate the quantities in the Geom struct in the MKSBHAC metric.
     */
    GRRT_SHARED_FUNCTION Geom make_geom(
      scalar a,
      scalar hslope,
      scalar x1, // position
      scalar x2  // wavevector
    ) {
      Geom G;
      G.a = a;
      G.a2 = a * a;
      G.r = std::exp(x1);
      G.r2 = G.r * G.r;
      scalar s2x, c2x;
      sincos_t(scalar(2) * x2, &s2x, &c2x);
      const scalar th = x2 + scalar(0.5) * hslope * s2x;
      G.tp = scalar(1) + hslope * c2x;
      G.tp2 = G.tp * G.tp;
      G.tpp = -scalar(2) * hslope * s2x;
      sincos_t(th, &G.s, &G.c);
      G.s2 = G.s * G.s;
      G.s2th = scalar(2) * G.s * G.c;
      G.cot = G.c / G.s;
      G.Sig = G.r2 + G.a2 * G.c * G.c;
      G.Sb = scalar(1) / G.Sig;
      G.Delta = G.r2 - scalar(2) * G.r + G.a2;
      G.A = G.Sig * G.Delta + scalar(2) * G.r * (G.r2 + G.a2);
      G.B = (scalar(2) * G.r2 - G.Sig) * G.Sb * G.Sb * G.Sb;
      G.C = G.r * G.Sb - G.a2 * G.B * G.s2;
      return G;
    }

    /**
     * @brief Symmetric bilinear form out[i] = -Gamma^i_jk U^j V^k.
     */
    GRRT_SHARED_FUNCTION void accel(
      const Geom* geom_in,
      const scalar U[4],
      const scalar V[4],
      scalar out[4]
    ) {
      const Geom G = *geom_in;
      const scalar h = scalar(0.5), two = scalar(2), four = scalar(4);
      const scalar P00 = U[0] * V[0];
      const scalar P01 = h * (U[0] * V[1] + U[1] * V[0]);
      const scalar P02 = h * (U[0] * V[2] + U[2] * V[0]);
      const scalar P03 = h * (U[0] * V[3] + U[3] * V[0]);
      const scalar P11 = U[1] * V[1];
      const scalar P12 = h * (U[1] * V[2] + U[2] * V[1]);
      const scalar P13 = h * (U[1] * V[3] + U[3] * V[1]);
      const scalar P22 = U[2] * V[2];
      const scalar P23 = h * (U[2] * V[3] + U[3] * V[2]);
      const scalar P33 = U[3] * V[3];

      // dr/dx1 = r for r = exp(x1), hence rp = r
      const scalar a = G.a, a2 = G.a2, r = G.r, r2 = G.r2, rp = G.r, rp2 = G.r2;
      const scalar tp = G.tp, tp2 = G.tp2, s = G.s, c = G.c, s2 = G.s2, s2th = G.s2th;
      const scalar Sig = G.Sig, Sb = G.Sb, Sb2 = Sb * Sb, Sb3 = Sb2 * Sb;
      const scalar Del = G.Delta, A = G.A, B = G.B, C = G.C, cot = G.cot;

      // t
      const scalar K = a2 * r * s2th * Sb2 * tp;
      out[0] = -two * (-r2 * Sb * tp2 * P22 - r * C * s2 * P33
                       + K * (-P02 - rp * P12 + a * s2 * P23)
                       + B * ((Sig + two * r) * rp * P01 + r * P00 + (Sig + r) * rp2 * P11
                              - two * a * r * s2 * P03 - a * s2 * (Sig + two * r) * rp * P13));

      // x1
      out[1] = -((scalar(1) - (two * Sig - Del) * B * rp) * P11
                 + (Del / rp) * (B * P00 - two * a * B * s2 * P03 - r * Sb * tp2 * P22 - C * s2 * P33)
                 + two * ((Del - Sig) * B * P01 - a2 * s * c * Sb * tp * P12
                          + a * (r * Sb + (Sig - Del) * B) * s2 * P13));

      // x2
      out[2] = -((G.tpp / tp - a2 * s * c * Sb * tp) * P22
                 + Sb3 * ((-a2 * r * s2th) * (P00 + rp2 * P11)
                          - s * c * (Del * Sig * Sig + two * r * (r2 + a2) * (r2 + a2)) * P33
                          - two * a2 * r * s2th * rp * P01
                          + two * a * s * c * (A + Sig * (Sig - Del)) * rp * P13
                          + two * a * r * (r2 + a2) * s2th * P03) / tp
                 + two * r * Sb * rp * P12);

      // x3
      out[3] = -(a * (B * P00 + B * rp2 * P11 - r * Sb * tp2 * P22 - C * s2 * P33)
                 + two * a * B * rp * P01
                 - four * a * r * cot * Sb2 * tp * P02
                 - two * a * (Sig + two * r) * cot * Sb2 * rp * tp * P12
                 - two * a2 * B * s2 * P03
                 + two * (cot + a2 * r * s2th * Sb2) * tp * P23
                 + two * C * rp * P13);
    }
  }
}