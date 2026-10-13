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

#pragma once

#include "core/geometry/headers.hpp"
#include "core/transfer/headers.hpp"

// Maximum number of observing frequencies carried in one pass.
// For each frequency, the cost is roughly 0.3 KB (FP32) to 0.6 KB (FP64)
#ifndef PHX_MAX_FREQS
  #define PHX_MAX_FREQS 32
#endif

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
    struct RunStats {
      long steps, nplasma;
      long bad;  // segments whose operator came out non-finite
    };

    /**
     * @brief Per-pixel state that survives from one sample to the next
     */
    struct PixelState {
      core::operators::Accumulator acc[PHX_MAX_FREQS];
      scalar tau_sum[PHX_MAX_FREQS];
      scalar tauF_sum[PHX_MAX_FREQS];  // diagnostic optical depths
    };

    GRRT_SHARED_FORCEINLINE void process_sample(
      const core::bhac::Model* model,
      PixelState* ps,
      const core::geodesic::Sample* sm,
      int right_leaf
    ) {
      
    }

    /**
     * @brief Renders one pixel
     *
     * @param model  Model
     * @param alpha  Horizontal impact parameter [Rg]
     * @param beta   Vertical impact parameter [Rg]
     * @param out    Output, nfreq entries
     * @return       Ray statistics
     */
    GRRT_SHARED_FUNCTION RunStats render_pixel(
      const core::bhac::Model* model,
      double alpha,
      double beta,
      FreqOut* out
    ) {
      const int nfreqs = model->uparams.nfreqs;

      RunStats stats;
      stats.steps = stats.nplasma = stats.bad = 0;
      for (int q = 0; q < nfreqs; q++) {
        for (int c = 0; c < 4; c++) {
          out[q].IQUV[c] = 0.;
        }
        out[q].tau = out[q].tauF = 0.;
      }

      // ray parameters (x_u and k_u)
      double photon_u[8];
      // Polarization basis
      double F1_u[4];
      // These are only used to populate the following two functions
      core::camera::initialize_photon(model, alpha, beta, photon_u);
      core::tetrad::camera_basis(model, photon_u, photon_u + 4, F1_u);

      PixelState ps;
      FOR_FREQ(f, nfreqs) {
        core::operators::init_accumulator(&ps.acc[f]);
        ps.tau_sum[f] = ps.tauF_sum[f] = scalar(0);
      }

      core::geodesic::Geodesic geo;
      core::geodesic::init(&geo, photon_u, F1_u);

      int steps = 0;
      int terminate = 0;
      scalar dl_prev = scalar(0);
      core::geodesic::Sample smp_pending;
      smp_pending.i = -1;

      while (true) {
        core::geodesic::Sample smp_curr;
        core::geodesic::view(&geo, &smp_curr);
        smp_curr.i = steps;
        smp_curr.dl = dl_prev;
        const scalar r = std::exp(smp_curr.x[1]);
        const int alive = r < model->r_out && r > model->cut_in && steps < model->uparams.max_steps && !terminate;
        if (smp_pending.i >= 1) {
          process_sample(model, &ps, &smp_curr, alive ? core::sampler::leaf_of(model, smp_curr.x) : -1);
        }


        break;
      }

      return stats;
    }
  }
}