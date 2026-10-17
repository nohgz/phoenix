/**
 * @file main.cpp
 * @author C. Noga (github.com/nohgz)
 * @brief Runs the show!
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "host/formatting/headers.hpp"
#include "host/io/headers.hpp"
#include "core/model/headers.hpp"
#include "host/render/headers.hpp"
#include "core/render/headers.hpp"

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads(void) { return 1; }
#endif

#if PHX_GPU_ENABLED
#include "device/headers.cuh"
#endif

namespace fp = host::fancy_printing;

static double now_s(void) {
  using clk = std::chrono::steady_clock;
  return std::chrono::duration<double>(clk::now().time_since_epoch()).count();
}


/**
 * @brief Writes the image (HDF5 and/or text) and the spectrum.
 */
static void write_output(
  const core::bhac::Model* model,
  const host::config::Params& params,
  const std::vector<core::pixel::FreqOut>& res
) {
  const int nfreqs = model->uparams.nfreqs;
  bool want_text = (params.image == "text" || params.image == "both");
  bool want_h5 = (params.image == "hdf5" || params.image == "both");

  if (want_h5) {
    const std::string hp = params.outdir + "/img_data.h5";
    if (host::output::write_image_hdf5(hp, model, params.block_size, res)) {
      fp::okmsg("Wrote %s\n", hp.c_str());
    } else {
      fp::warnmsg("Failed to write \"%s\" writing text instead...", hp.c_str());
      want_text = true;
    }
  }
  if (want_text) {
    for (int f = 0; f < nfreqs; f++) {
      const std::string ip = params.outdir + "/image" + (nfreqs == 1 ? "" : "_" + std::to_string(f)) + ".dat";
      host::output::write_image(ip, model, res, f);
      fp::okmsg("Wrote %s\n", ip.c_str());
    }
  }
  const std::string sp = params.outdir + "/spectrum" + ".dat";
  host::output::write_spectrum(sp, model, res);
  fp::okmsg("Wrote %s", sp.c_str());
  fp::hline();
  for (int f = 0; f < nfreqs; f++) {
    double S[4];
    host::output::compute_spec(model, res, f, S);
    fp::infomsg("Flux at %.4e Hz (Jy): \033[1;38;5;160mI\033[0m %+.9e \033[1;38;5;77mQ\033[0m %+.9e \033[1;38;5;27mU\033[0m %+.9e \033[1;38;5;11mV\033[0m %+.9e\n", model->uparams.frequencies[f], S[0], S[1], S[2], S[3]);
  }
}

int main(int argc, char *argv[]) {
  const double t_start = now_s();

  std::string cfgpath = "input/config.toml";

  // -- parse cli arguments, which can override the config.toml --
  int first_arg_idx = 1;
  // if the first arg is the config.toml filepath
  if (argc > 1 && strchr(argv[1], '=') == NULL) {
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
      fp::infomsg("usage: %s [input/config.toml] [section.key=value ...]\n", argv[0]);
      return 0;
    }
    cfgpath = argv[1];
    first_arg_idx = 2;
  }
  std::vector<std::string> cli_overrides(argv + first_arg_idx, argv + argc);

#if PHX_GPU_ENABLED
  device::warmup_start();  // CUDA context on a background thread while the config and the dump are read
#endif

  core::bhac::Model model = {};
  host::config::Params conf_params;

  try {
    host::config::load_config(cfgpath, cli_overrides, conf_params, &model.uparams);
  } catch(const std::exception& e) {
    fp::failmsg("%s", e.what());
    return EXIT_FAILURE;
  }

  if (conf_params.show_header) {
    fp::header();
    fp::hline();
  }
  fp::infomsg("Input Filepaths\n");
  fp::okmsg("Load config file      \"%s\"", cfgpath.c_str());

#if !PHX_GPU_ENABLED
  if (conf_params.gpu_decode) {
    fp::warnmsg("No GPU available, falling back to host decoding!");
    conf_params.gpu_decode = false;
  }
#endif

  try {
    fp::blankmsg("Load grid information \"%s\"", conf_params.grid_file.c_str());
    host::bhac_io::load_bhac_grid(&model.gparams, conf_params.grid_file);
    fp::return_ok();
  } catch (const std::exception& e) {
    fp::failmsg("%s", e.what());
    return EXIT_FAILURE;
  }

  try {
    fp::blankmsg("Load GRMHD data       \"%s\"", conf_params.dump_file.c_str());
    const double t_model = now_s();
    host::bhac_io::load(&model, conf_params.dump_file, conf_params.gpu_decode);
    fp::return_ok();
    fp::timemsg(" ⇥ Loaded in %.2f s%s", now_s() - t_model, conf_params.gpu_decode ? "   \033[1;33m(decoding deferred to GPU)\033[0m" : "");
  } catch (const std::exception& e) {
    fp::failmsg("%s", e.what());
    return EXIT_FAILURE;
  }

  const core::bhac::UserParams& up = model.uparams;
  fp::hline();
  fp::infomsg("Precision: %s\n", accuracy_name());
  fp::infomsg("Geometry:\n     ⇥ a       = %g\n     ⇥ hslope  = %g\n     ⇥ Horizon = %g Rg\n",double(model.a), double(model.hslope), double(model.cut_in));
  fp::infomsg("Black Hole & Plasma:\n     ⇥ Mass    = %.4e M_sun\n     ⇥ Dist    = %.4e kpc\n     ⇥ M_unit  = %.2e grams\n     ⇥ r_low   = %g\n     ⇥ r_high  = %g\n     ⇥ beta_c  = %g\n", up.mass, up.distance, up.m_unit, up.r_low, up.r_high, up.beta_c);
  fp::infomsg("Camera:\n     ⇥ Inclination = %g deg\n     ⇥ Azimuth     = %g deg\n     ⇥ Pos. Angle  = %g deg\n     ⇥ Resolution  = %d x %d px\n     ⇥ FOV         = %g x %g Rg\n     ⇥ r_cam       = %g Rg\n     ⇥ Frequencies = %d", up.inclination, up.phi_cam, up.pa, up.img_width, up.img_height, up.fov_x, up.fov_y, up.r_cam, up.nfreqs);

  std::vector<core::pixel::FreqOut> res((size_t)up.img_width * up.img_height * up.nfreqs);

  const std::vector<int> jobs = host::jobs::schedule(&up, PHX_GPU_ENABLED);
  const int njobs = int(jobs.size());
  std::atomic<int> ndone(0);
  core::pixel::RunStats tot = {0, 0, 0};

  fp::hline();
  fprintf(stderr, "\n");
#if PHX_GPU_ENABLED
  device::Info gi;
#endif
  const double t_ray = now_s();
  try {
    {
      auto bar = barkeep::ProgressBar(&ndone, {
        .out = &std::cerr,
        .total = njobs,
        .message = "Raymarching",
        .speed = 0.1,
        .speed_unit = "px/s",
        .style = barkeep::Rich,
        .no_tty = !isatty(fileno(stderr))
      });
      try {
#if PHX_GPU_ENABLED
        device::render(&model, conf_params.dump_file.c_str(), jobs.data(), njobs, res.data(), &gi, device::count_progress, &ndone);
#else
        tot = host::jobs::process(&model, jobs, res.data(), &ndone);
#endif
      } catch (...) {
        bar->done();
        throw;
      }
      bar->done();
    }
#if PHX_GPU_ENABLED
    tot.steps = gi.steps, tot.nplasma = gi.nplasma, tot.bad = gi.bad;
#endif
  } catch (const std::exception& e) {
    fp::failmsg("%s", e.what());
    host::bhac_io::free_model(&model);
    return EXIT_FAILURE;
  }
  fp::timemsg("Raymarch: %.2f s\n", now_s() - t_ray);
#if PHX_GPU_ENABLED
  if (gi.decode_s > 0) {
    fp::infomsg(" ⇥ Model upload %.2f s, GPU decode: %.2f s, kernel %.2f s\n", gi.upload_s, gi.decode_s, gi.kernel_s);
  } else {
    fp::infomsg(" ⇥ Model upload %.2f s, kernel %.2f s\n", gi.upload_s, gi.kernel_s);
  }
#endif
  fp::hline();

  long badpix = 0;
  for (size_t p = 0; p < (size_t)up.img_width * up.img_height; p++) {
    int nonfinite = 0;
    for (int f = 0; f < up.nfreqs; f++) {
      for (int q = 0; q < 4; q++) nonfinite |= !std::isfinite(res[p * up.nfreqs + f].IQUV[q]);
    }
    badpix += nonfinite;
  }
  fp::infomsg("%.0f geodesic samples/pixel, %.0f plasma segments/pixel (mean)\n", (double)tot.steps / njobs, (double)tot.nplasma / njobs);
  fp::infomsg("Non-finite segment operators: %ld, pixels with non-finite IQUV: %ld\n", tot.bad, badpix);
#if PHX_GPU_ENABLED
  fp::infomsg("Device: %s, %d SMs\n", gi.device, gi.sms);
  fp::infomsg("Launch: %d blocks x %d threads, %d registers/thread, %d B local memory/thread", gi.blocks, gi.block, gi.regs, gi.local_bytes);
#endif

  fp::hline();
  const double t_out = now_s();
  write_output(&model, conf_params, res);
  fp::timemsg("Output: %.2f s\n", now_s() - t_out);

  host::bhac_io::free_model(&model);
  fp::timemsg("Execution: %.2f s\n", now_s() - t_start);
  return EXIT_SUCCESS;
}
