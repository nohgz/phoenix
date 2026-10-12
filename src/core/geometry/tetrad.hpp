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
     * @brief levi-civita permutation symbol with eps[0,1,2,3] = 1
     */
    GRRT_SHARED_FUNCTION int lc_eps(int i, int j, int k, int l) {
      if (i == j || i == k || i == l || j == k || j == l || k == l) return 0;
      return ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l)) / 12;
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

      scalar perp[4];
      for (int i = 0; i < 4; i++) {
        perp[i] = scalar(0);
        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
              const int e = lc_eps(i, j, k, l);
              if (e) perp[i] += (scalar(-1) / sq) * scalar(e) * U_d[j] * k_d[k] * b_d[l] / (omega * Nc);
            }
          }
        }
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
  }
}