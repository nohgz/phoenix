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

#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "host/formatting/headers.hpp"
#include "host/io/headers.hpp"
#include "core/geometry/headers.hpp" // for eps4 test
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

int main(int argc, char *argv[]) {

  std::string cfgpath = "input/config.toml";

  // -- parse cli arguments, which can override the config.toml --
  int first_arg_idx = 1;
  // if the first arg is the config.toml filepath
  if (argc > 1 && strchr(argv[1], '=') == NULL) {
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
      host::fancy_printing::infomsg("usage: %s [input/config.toml] [section.key=value ...]\n", argv[0]);
      return 0;
    }
    cfgpath = argv[1];
    first_arg_idx = 2;
  }
  std::vector<std::string> cli_overrides(argv + first_arg_idx, argv + argc);

  core::bhac::Model model = {};
  host::config::Params conf_params;

  try {
    host::config::load_config(cfgpath, cli_overrides, conf_params, &model.uparams);
  } catch(const std::exception& e) {
    host::fancy_printing::failmsg(e.what());
    return EXIT_FAILURE;
  }

  if (conf_params.show_header) {
    host::fancy_printing::header();
    host::fancy_printing::hline();
  }
  host::fancy_printing::infomsg("Input Filepaths\n");
  host::fancy_printing::okmsg("Load config file      \"%s\"", cfgpath.c_str());

#if !PHX_GPU_ENABLED
  if (conf_params.gpu_decode) {
    host::fancy_printing::warnmsg("No GPU available, falling back to host decoding!");
    conf_params.gpu_decode = false;
  }
#endif

  try {
    host::fancy_printing::blankmsg("Load grid information \"%s\"", conf_params.grid_file.c_str());
    host::bhac_io::load_bhac_grid(&model.gparams, conf_params.grid_file);
    host::fancy_printing::return_ok();
  } catch (std::runtime_error e) {
    host::fancy_printing::failmsg(e.what());
    return EXIT_FAILURE;
  }

  try {
    host::fancy_printing::blankmsg("Load GRMHD data       \"%s\"", conf_params.dump_file.c_str());
    host::bhac_io::load(&model, conf_params.dump_file, conf_params.gpu_decode);
    host::fancy_printing::return_ok();
  } catch (std::runtime_error e) {
    host::fancy_printing::failmsg(e.what());
    return EXIT_FAILURE;
  }

  host::fancy_printing::hline();
  host::fancy_printing::infomsg("uhhh");
  host::fancy_printing::blankmsg("Create output directory \"%s\"", conf_params.outdir.c_str());
  mkdir(conf_params.outdir.c_str(), 0755);
  host::fancy_printing::return_ok();

  std::vector<core::pixel::FreqOut> res(
    (size_t)model.uparams.img_width * model.uparams.img_height * model.uparams.frequencies.size()
  );

  const std::vector<int> jobs = host::jobs::schedule(&model.uparams, PHX_GPU_ENABLED);
  const int njobs = (int)jobs.size();
  std::atomic<int> ndone(0);
  core::pixel::RunStats tot;

  fprintf(stderr, "\n");
  #if PHX_GPU_ENABLED
  #else
    host::fancy_printing::infomsg("Workload: %d pixels, %d thread(s), %0.1f px/thread\n\n", njobs, omp_get_max_threads(), float(njobs)/omp_get_max_threads());
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
      tot = host::jobs::process(&model, jobs, res.data(), &ndone);
      bar->done();
    }
  #endif

  return EXIT_SUCCESS;
}