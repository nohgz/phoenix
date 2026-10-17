/**
 * @file kernel.cuh
 * @author C. Noga (github.com/nohgz)
 * @brief The host-callable interface of the CUDA raymarcher
 *
 * GNU GPL v3.0
 */

#pragma once

#include <stddef.h>

#include "core/model/headers.hpp"
#include "core/render/headers.hpp"

namespace device {
  using core::bhac::Model;
  using core::pixel::FreqOut;

  /**
   * @brief What a render() call reports.
   */
  struct Info {
    char device[256];
    int sms, blocks, block;     // launch geometry (persistent grid = blocks, threads per block = block)
    int regs, local_bytes;      // registers / thread, local memory / thread (bytes)
    double upload_s, kernel_s;  // time for model upload, kernel time
    double decode_s;            // file read + GPU decode (grmhd.gpu_decode)
    long steps, nplasma, bad;   // sums over all rays for diagnostics
  };

  using Progress = void (*)(int jobs_done_in_chunk, void* user);

  /**
   * @brief Adds n to the std::atomic<int> that user points to. Used in the progress bar.
   */
  void count_progress(int n, void* user);

  /**
   * @brief CUDA context creation costs about 0.1-0.2 s on the first runtime call.
   *        This creates it on a background thread
   */
  void warmup_start();

  /**
   * @brief Renders a list of pixels on the GPU.
   *
   * @param model    Model on the host (with prim == NULL the GPU reads and decodes the dump itself)
   * @param dump     The dump file (only read when m->prim is NULL)
   * @param jobs     Pixel indices i * img_h + j, ideally Morton ordered (consecutive 32 = one warp's tile)
   * @param njobs    Number of jobs
   * @param res      img_w * img_h * nfreq FreqOut, res[(i*img_h + j)*nfreq + f]
   * @param info     Output
   * @param progress Called after each chunk of the job list (may be NULL)
   * @param user     Passed to progress
   */
  void render(
    const Model* model,
    const char* dump,
    const int* jobs,
    int njobs,
    FreqOut* res,
    Info* info,
    Progress progress,
    void* user
  );
}
