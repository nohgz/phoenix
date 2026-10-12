/**
 * @file pixel.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief
 * @version 1.0
 * @date 2026-10-10
 *
 * GNU GPL v3.0
 *
 */

#include "core/geometry/headers.hpp"

namespace core {
  namespace pixel {
    /**
     * @brief Result for one pixel and one frequency.
     */
    struct FreqOut {
      double IQUV[4];
      double tau, tauF;
    };

    /**
     * @brief Ray statistics for the run.
     */
    struct PixelStats {
      long steps, nplasma;
      long bad;  // segments whose operator came out non-finite
    };

    /**
     * @brief Renders one pixel
     *
     * @param model  Model
     * @param alpha  Horizontal impact parameter [Rg]
     * @param beta   Vertical impact parameter [Rg]
     * @param out    Output, nfreq entries
     * @return       Ray statistics
     */
    GRRT_SHARED_FUNCTION PixelStats render_pixel(
      const core::bhac::Model* model,
      double alpha,
      double beta,
      FreqOut* out
    ) {
      const int nf = model->uparams.nfreqs;

      PixelStats stats;
      stats.steps = stats.nplasma = stats.bad = 0;
      for (int q = 0; q < nf; q++) {
        for (int c = 0; c < 4; c++) {
          out[q].IQUV[c] = 0.;
        }
        out[q].tau = out[q].tauF = 0.;
      }

      // Current ray parameters (x_u and k_u)
      double photon_u[8];
      // Polarization basis
      double F1_u[4];

      core::camera::initialize_photon(model, alpha, beta, photon_u);
      core::tetrad::camera_basis(model, photon_u, photon_u + 4, F1_u);

      return stats;
    }
  }
}