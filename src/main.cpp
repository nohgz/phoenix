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

#include "host/formatting/headers.hpp"
#include "host/io/headers.hpp"
#include "core/geometry/headers.hpp" // for eps4 test
#include "core/model/headers.hpp"

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
    host::fancy_printing::infomsg("No GPU available, falling back to host decoding!");
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

  #if PHX_GPU_ENABLED
    fprintf(stderr, "\nYep, we have a GPU!\n");
    device::launch_hello();
  #else
    fprintf(stderr, "\nNah, we don't have a GPU!\n");
    fprintf(stderr, "[CPU] eps4(1,2,3,4) returned result %d\n", core::geometry::eps4(1, 2, 3, 4));
  #endif

  return EXIT_SUCCESS;
}