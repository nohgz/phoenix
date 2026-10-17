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
#include "core/model/headers.hpp"
#include "core/transfer/headers.hpp"

// Maximum number of observing frequencies carried in one pass.
// For each frequency, the cost is roughly 0.3 KB (FP32) to 0.6 KB (FP64)

namespace core {
  namespace pixel {
    /**
     * @brief Result for one pixel and one frequency.
     */
    struct FreqOut {
      double IQUV[4]; // erg/s/cm^2/Hz/sr (Jy / 1e23)
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
      int nplasma, bad;                // more diagnostics
    };

    /**
     * @brief The meat and potatoes of everything
     */
    GRRT_SHARED_FORCEINLINE void process_sample(
      const core::bhac::Model* model,
      PixelState* ps,
      const core::geodesic::Sample* sm,
      int right_leaf
    ) {
      const int nfreqs = model->uparams.nfreqs;

      scalar X_u[4], k_u[4], F1_u[4];
      for (int q = 0; q < 4; q++) {
        X_u[q] = sm->x[q];
        k_u[q] = sm->k[q];
        F1_u[q] = sm->F1[q];
      }
      const scalar r = std::exp(X_u[1]);
      scalar g_dd[4][4];
      core::metric::metric_dd<scalar>(model->a, model->hslope, X_u, g_dd);
      if (std::fabs(core::metric::inner_product(g_dd, k_u, k_u)) > scalar(1e-6) && r > scalar(2)) {
        core::geodesic::normalize_null(g_dd, k_u);
      }

      core::sampler::Fluid fl = {};
      int cache = right_leaf;
      const int in_volume = core::sampler::get_fluid_params(model, X_u, &fl, &cache);

      // make the geometry of the sample in the plasma frame
      int emit = 0;
      scalar ratio = scalar(0), c2 = scalar(0), s2 = scalar(0), pitch = scalar(0);
      if (in_volume && r < model->uparams.rt_outer) {
        scalar smp_g_dd[4][4], e_u[4][4];
        core::metric::metric_dd<scalar>(model->a, model->hslope, fl.X_u, smp_g_dd);  // the position the sampler used
        pitch = core::tetrad::pitch_angle(smp_g_dd, k_u, fl.B_u, fl.U_u);
        if (!(pitch == scalar(0) || pitch == scalar(M_PI))) {
          core::tetrad::create_observer_tetrad(model, fl.X_u, k_u, fl.U_u, fl.B_u, e_u);
          core::tetrad::chi_from_F1(smp_g_dd, e_u, F1_u, &c2, &s2);  // only F1 is transported
          ratio = -core::metric::inner_product(smp_g_dd, k_u, fl.U_u);     // nu_plasma / nu_observed
          emit = 1;
        }
      }
      if (!emit) return;

      // do this segment per frequency. coefficients and operator are swept in
      // one pass, so we only store the accumulator matrices
      FOR_FREQ(f, nfreqs) {
        const scalar nu_p = model->freq[f] * ratio;
        const core::coeffs::PolCoeffs raw = core::coeffs::thermal_raw(fl.theta_e, fl.n_e, nu_p, fl.B, pitch);
        const core::coeffs::PolCoeffs inv = core::coeffs::invariant(raw, nu_p, model->jscale[f]);
        const core::coeffs::PolCoeffs seg = core::coeffs::rotate_to_camera(&inv, c2, s2);

        const scalar h = sm->dl * model->hconst[f];  // dl * (m_e c^2 / h nu) * r_g
        scalar K[4][4], j[4];
        core::operators::build_K(&seg, K);
        core::operators::emis_vec(&seg, j);
        const core::operators::Accumulator op = core::operators::trapezoid(K, j, h);

        int finite = 1;
        for (int a = 0; a < 4; a++) {
          finite = finite && std::isfinite(op.S[a]);
          for (int b = 0; b < 4; b++) finite = finite && std::isfinite(op.A[a][b]);
        }
        // Drop the segment instead of poisoning the pixel with NaNs
        if (!finite) {
          ps->bad++;
          continue;
        }
        core::operators::accumulate(&ps->acc[f], &op);

        ps->tau_sum[f] += seg.aI * h;
        ps->tauF_sum[f] += std::fabs(seg.rV * h);  // diagnostics
      }
      ps->nplasma++;
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

      ps.nplasma = ps.bad = 0;

      int steps = 0, turns = 0, terminate = 0;
      scalar thdot_prev = scalar(0), dl_prev = scalar(0);
      const int beta_neg = beta < 0., beta_pos = beta > 0.;
      core::geodesic::Sample smp_pending = {};  // waits for its far neighbour
      smp_pending.i = -1;

      // One call site for process_sample. Each pass visits the current sample, then (if the ray goes on) steps.
      // A sample is processed one pass late, when the sampler's cached block is the one of its far neighbour;
      // the final pass (ray ended) flushes it with no cached block.
      for (;;) {
        core::geodesic::Sample smp_curr;
        core::geodesic::view(&geo, &smp_curr);
        smp_curr.i = steps;     // state recorded BEFORE the next step
        smp_curr.dl = dl_prev;  // the step that arrived here
        const scalar r = std::exp(smp_curr.x[1]);
        const int alive = r < model->r_out && r > model->cut_in && steps < model->uparams.max_steps && !terminate;
        if (smp_pending.i >= 1) {
          process_sample(model, &ps, &smp_pending, alive ? core::sampler::leaf_of(model, smp_curr.x) : -1);
        }
        smp_pending = smp_curr;
        if (!alive) break;

        // turning points of theta (image order)
        if (thdot_prev * smp_curr.k[2] < scalar(0) && steps > 2) turns++;
        thdot_prev = smp_curr.k[2];
        if ((beta_neg && turns > model->uparams.max_order) || (beta_pos && turns > (model->uparams.max_order + 1))) {
          terminate = 1;
        }

        // grid-aware step size, then RK4
        const scalar* dxc = core::sampler::grid_dxc(model, smp_curr.x);
        scalar X_u[4], k_u[4];
        for (int i = 0; i < 4; i++) {
          X_u[i] = smp_curr.x[i];
          k_u[i] = smp_curr.k[i];
        }
        X_u[2] = core::geodesic::pole_x2(&geo);  // distance to the nearer pole when folding (FP32), else x2 itself
        const scalar dl = core::geodesic::step_size(model->uparams.stepsize, X_u, k_u, dxc);
        core::geodesic::step(model->a, model->hslope, &geo, dl);
        dl_prev = std::fabs(dl);
        steps++;
      }

      stats.steps = steps;
      stats.nplasma = ps.nplasma;
      stats.bad = ps.bad;
      if (ps.nplasma > 0) {
        FOR_FREQ(f, nfreqs) {
          out[f].tau = (double)ps.tau_sum[f];
          out[f].tauF = (double)ps.tauF_sum[f];
          for (int q = 0; q < 4; q++) out[f].IQUV[q] = (double)ps.acc[f].S[q];  // already scaled by nu^3
        }
      }
      return stats;
    }
  }
}