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

#if PHX_GPU_ENABLED
#include "device/headers.cuh"
#endif

int main(int argc, char *argv[]) {
  
  std::string cfgpath = "input/config.toml";

  // -- parse cli arguments, which can override the config.toml --
  int first_arg_idx = 1;
  if (argc > 1 && strchr(argv[1], '=') == NULL) {
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
      host::fancy_printing::infomsg("usage: %s [input/config.toml] [section.key=value ...]\n", argv[0]);
      return 0;
    }
    cfgpath = argv[1];
    first_arg_idx = 2;
  }
  std::vector<std::string> cli_overrides(argv + first_arg_idx, argv + argc);

  host::config::Config conf;

  try {
    host::config::load_config(cfgpath, cli_overrides, conf);
  } catch(const std::exception& e) {
    host::fancy_printing::failmsg(e.what());
    return EXIT_FAILURE;
  }
  
  if (conf.show_header) host::fancy_printing::header();

  #if PHX_GPU_ENABLED
    fprintf(stderr, "Yep, we have a GPU!\n");
    device::launch_hello();
  #else 
    fprintf(stderr, "Nah, we don't have a GPU!\n"); 
    fprintf(stderr, "[CPU] eps4(1,2,3,4) returned result %d\n", core::geometry::eps4(1, 2, 3, 4));
  #endif

  return EXIT_SUCCESS;
}