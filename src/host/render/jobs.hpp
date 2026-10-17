#pragma once

#include <algorithm>
#include <atomic>
#include <utility>
#include <vector>

#include "host/io/headers.hpp"
#include "core/render/headers.hpp"

namespace host {
  namespace jobs {
    /**
     * @brief Give every pixel an index i * img_height + j. If on the GPU, the
     * list is in Morton (Z-order), so that consecutive 32 (1 CUDA warp!) jobs
     * are an 8x4 / 4x8 tile of pixels.
     */
    static inline std::vector<int> schedule(
      const core::bhac::UserParams* uparams,
      bool morton
    ) {
      std::vector<int> jobs;
      for (int i = 0; i < uparams->img_width; i++) {
        for (int j = 0; j < uparams->img_height; j++) {
          jobs.push_back(i * uparams->img_height + j);
        }
      }
      if (morton) {
        std::vector<std::pair<unsigned long long, int>> keyed(jobs.size());
        for (size_t n = 0; n < jobs.size(); n++) {
          const int i = jobs[n] / uparams->img_height;
          const int j = jobs[n] % uparams->img_height;
          unsigned long long z = 0;
          for (int b = 0; b < 16; b++) {
            z |= ((unsigned long long)((i >> b) & 1) << (2 * b)) | ((unsigned long long)((j >> b) & 1) << (2 * b + 1));
          }
          keyed[n] = std::make_pair(z, jobs[n]);
        }
        std::sort(keyed.begin(), keyed.end());
        for (size_t n = 0; n < jobs.size(); n++) {
          jobs[n] = keyed[n].second;
        }
      }
      return jobs;
    }

    /**
     * @brief Renders the pixels in jobs on all available OpenMP threads
     */
    inline core::pixel::RunStats process(
      const core::bhac::Model* model,
      const std::vector<int>& jobs,
      core::pixel::FreqOut* result,
      std::atomic<int>* ndone
    ) {
      long steps = 0, segs = 0, bad = 0;
      const int njobs = (int)jobs.size();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : steps, segs, bad)
      for (int n = 0; n < njobs; n++) {
        const int idx = jobs[n];
        const int i = idx / model->uparams.img_height;
        const int j = idx % model->uparams.img_height;

        double a0, b0;
        core::camera::pixel_impact(
          i,
          j,
          model->uparams.img_width,
          model->uparams.img_height,
          model->uparams.fov_x,
          model->uparams.fov_y,
          &a0, &b0
        );
        // position angle: the ray of the unrotated sky that lands on this pixel (cos 0 = 1, sin 0 = 0: none)
        const double alpha = a0 * model->cos_pa + b0 * model->sin_pa;
        const double beta = -a0 * model->sin_pa + b0 * model->cos_pa;

        const core::pixel::RunStats ps = core::pixel::render_pixel(
          model,
          alpha,
          beta,
          &result[(size_t)idx * model->uparams.nfreqs]
        );

        steps += ps.steps;
        segs += ps.nplasma;
        bad += ps.bad;
        ++*ndone;
      }
      core::pixel::RunStats total = {steps, segs, bad};
      return total;
    }
  }
}