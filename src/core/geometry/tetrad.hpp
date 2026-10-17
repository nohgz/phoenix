/**
 * @file tetrad.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief
 * @version 1.0
 * @date 2026-10-11
 *
 * GNU GPL v3.0
 *
 */
#pragma once

#include "core/numerics/headers.hpp"
#include "core/model/headers.hpp"
#include "core/geometry/metric.hpp"

namespace core {
  namespace tetrad {

  /**
   * @brief U^mu of the camera observer
   * @param model
   * @param X0_u  Camera position
   * @param k0_u  Camera wavevector
   * @param F1_u  Output polarization vector
   */
    GRRT_SHARED_FUNCTION void construct_U(
      const core::bhac::Model* model,
      const scalar X_u[4],
      scalar U_u[4]
    ) {
      scalar g_uu[4][4];
      core::metric::metric_uu<scalar>(model->a, model->hslope, X_u, g_uu);
      const scalar Ud0 = scalar(-1);
      const scalar B = -g_uu[0][3] * Ud0 / g_uu[3][3];
      const scalar C = -(scalar(1) + g_uu[0][0] * Ud0 * Ud0) / g_uu[3][3];
      scalar U_d[4] = {Ud0, scalar(0), scalar(0), B + std::sqrt(B * B + C)};
      core::metric::raise_index(g_uu, U_d, U_u);
    }

    /**
     * @brief RAPTOR create_observer_tetrad()
     *
     * @param m    Model
     * @param X_u  Position
     * @param k_u  Wavevector
     * @param U_u  Observer 4-velocity
     * @param b_u  Reference vector
     * @param e_u  Output tetrad e_u
     */
    GRRT_SHARED_FORCEINLINE void create_observer_tetrad(
      const core::bhac::Model* m,
      const scalar X_u[4],
      const scalar k_u[4],
      const scalar U_u[4],
      const scalar b_u[4],
      scalar e_u[4][4]
    ) {
      scalar g_dd[4][4];
      core::metric::metric_dd<scalar>(m->a, m->hslope, X_u, g_dd);

      const scalar omega = -core::metric::inner_product(g_dd, k_u, U_u);
      scalar eK[4];
      for (int i = 0; i < 4; i++) {
        eK[i] = k_u[i] / omega - U_u[i];
      }

      const scalar Beta = core::metric::inner_product(g_dd, U_u, b_u);
      const scalar Cc = core::metric::inner_product(g_dd, k_u, b_u) / omega - Beta;
      const scalar b2 = core::metric::inner_product(g_dd, b_u, b_u);
      const scalar Nc = std::sqrt(b2 + Beta * Beta - Cc * Cc);

      scalar para[4];
      for (int i = 0; i < 4; i++) {
        para[i] = (b_u[i] + Beta * U_u[i] - Cc * eK[i]) / Nc;
      }

      scalar U_d[4], k_d[4], b_d[4];
      core::metric::lower_index(g_dd, U_u, U_d);
      core::metric::lower_index(g_dd, k_u, k_d);
      core::metric::lower_index(g_dd, b_u, b_d);
      const scalar sq = std::sqrt(-core::metric::det4(g_dd));

      scalar perp[4] = {scalar(0), scalar(0), scalar(0), scalar(0)};
      const scalar s = scalar(-1) / (sq * omega * Nc);
      static const int idx[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
      GRRT_UNROLL for (int i = 0; i < 4; i++) {
        const int a = idx[i][0], c = idx[i][1], d = idx[i][2];
        const scalar det = U_d[a] * (k_d[c]*b_d[d] - k_d[d]*b_d[c])
                        - U_d[c] * (k_d[a]*b_d[d] - k_d[d]*b_d[a])
                        + U_d[d] * (k_d[a]*b_d[c] - k_d[c]*b_d[a]);
        perp[i] = ((i & 1) ? -s : s) * det;
      }

      for (int i = 0; i < 4; i++) {
        e_u[i][0] = U_u[i];
        e_u[i][1] = perp[i];
        e_u[i][2] = para[i];
        e_u[i][3] = eK[i];
      }
    }

    /**
     * @brief Camera polarization basis F1 (RAPTOR construct_f_obs_tetrad_u)
     *
     * @param model
     * @param X0_u  Camera position
     * @param k0_u  Camera wavevector
     * @param F1_u  Output polarization vector
     */
    GRRT_SHARED_FUNCTION void camera_basis(
      const core::bhac::Model* model,
      const double X0_u[4],
      const double k0_u[4],
      double F1_u[4]
    ) {
      scalar X_u[4], k_u[4], U_u[4], e_u[4][4];
      for (int i = 0; i < 4; i++) {
        X_u[i] = scalar(X0_u[i]);
        k_u[i] = scalar(k0_u[i]);
      }
      construct_U(model, X_u, U_u);
      const scalar cam_up[4] = {
        scalar(0),
        scalar(0),
        scalar(0),
        scalar(-1)
      };
      create_observer_tetrad(model, X_u, k_u, U_u, cam_up, e_u);
      for (int i = 0; i < 4; i++) {
        F1_u[i] = (double)e_u[i][1];
      }
    }

    /**
     * @brief The angle between k and B in the plasma frame
     */
    GRRT_SHARED_FORCEINLINE scalar pitch_angle(
      const scalar g_dd[4][4],
      const scalar k_u[4],
      const scalar B_u[4],
      const scalar U_u[4]
    ) {
      const scalar Bn = std::sqrt(std::fabs(core::metric::inner_product(g_dd, B_u, B_u)));
      if (Bn == scalar(0)) {
        return scalar(M_PI / 2.);
      }
      const scalar kk = std::fabs(core::metric::inner_product(g_dd, k_u, U_u));
      scalar mu = core::metric::inner_product(g_dd, k_u, B_u) / (kk * Bn);
      if (std::fabs(mu) > scalar(1)) {
        mu /= std::fabs(mu);
      }
      return std::acos(mu);
    }

    /**
     * @brief Rotation of the plasma -> camera Stokes map
     *
     * @param g_dd  Metric at the sample
     * @param e_u   Tetrad of the sample
     * @param F1_u  Transported camera polarization vector
     * @param c2    Output, cos 2psi
     * @param s2    Output, sin 2psi
     */
    GRRT_SHARED_FORCEINLINE void chi_from_F1(
      const scalar g_dd[4][4],
      const scalar e_u[4][4],
      const scalar F1_u[4],
      scalar* c2,
      scalar* s2
    ) {
      scalar ep[4], ea[4];
      for (int mu = 0; mu < 4; mu++) {
        ep[mu] = e_u[mu][1];
        ea[mu] = e_u[mu][2];
      }
      const scalar C = core::metric::inner_product(g_dd, ep, F1_u);
      const scalar S = core::metric::inner_product(g_dd, ea, F1_u);
      const scalar n = C * C + S * S;
      *c2 = (C * C - S * S) / n;
      *s2 = scalar(2) * C * S / n;
    }
  }
}