/**
 * @file config.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains the machinery that is used to read the TOML config file.
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/model/headers.hpp"

namespace host {
  namespace config {
    /**
    * @brief The code execution parameters passed in config.toml.
    *        See the comments in the config.toml file for descriptions of
    *        its members.
    */
    struct Params {
      // [visual]
      bool show_header;

      // [output]
      std::string outdir;
      // "hdf5" (RAPTOR format), "text", "both" or "none"
      std::string image = "hdf5";
      // pixels per block side in the HDF5 file; must divide the image size
      int block_size = 10;

      // [grmhd]
      std::string dump_file;
      std::string grid_file;
      bool gpu_decode;
    };


    /**
     * @brief Reads the observing frequencies and places them in the out vector
     *
     * @param kv    The key value pairs from the toml
     * @param used  Keys that were consumed
     * @param out   Output
     * @throw std::runtime_error
     */
    inline void read_frequencies(
      const std::map<std::string, std::string>& kv,
      std::set<std::string>& used,
      std::vector<double>& out
    ) {
      enum { SINGLE, LIST, LINSPACE, LOGSPACE, NFORMS };
      const char* forms[NFORMS] = {
        "camera.frequency",
        "camera.frequency_list",
        "camera.frequency_linspace",
        "camera.frequency_logspace"
      };

      int given = 0, form = SINGLE;
      for (int i = 0; i < NFORMS; i++) {
        if (kv.count(forms[i]) > 0) {
          given++;
          form = i;
        }
      }
      if (given != 1) {
        throw std::runtime_error("Too many (or no) frequency forms given!");
      }
      const std::string key = forms[form];
      const std::string& text = kv.at(key);
      used.insert(key);

      auto number = [&](const std::string& t) {
        char* end = nullptr;
        const double v = std::strtod(t.c_str(), &end);
        if (end == t.c_str() || *end != '\0') {
          throw std::runtime_error(key + ": not a number: " + t);
        }
        return v;
      };

      auto numbers = [&]() {
        const size_t a = text.find('['), b = text.rfind(']');
        if (a == std::string::npos || b == std::string::npos || b < a) {
          throw std::runtime_error(key + " must be a list on ONE line: [a, b, ...]");
        }
        std::vector<double> v;
        std::stringstream ss(text.substr(a + 1, b - a - 1));
        std::string tok;
        while (std::getline(ss, tok, ',')) {
          const size_t p = tok.find_first_not_of(" \t"), q = tok.find_last_not_of(" \t");
          if (p == std::string::npos) continue;  // empty (trailing comma)
          v.push_back(number(tok.substr(p, q - p + 1)));
        }
        return v;
      };

      out.clear();
      if (form == SINGLE) {
        out.push_back(number(text));
      } else {
        const std::vector<double> v = numbers();
        if (form == LIST) {
          if (v.empty()) {
            throw std::runtime_error(key + " is empty!");
          }
          out = v;
        } else {
          if (v.size() != 3) {
            throw std::runtime_error(key + " has " + std::to_string(v.size()) + "arguments! (Expected 3)");
          }
          const double lo = v[0], hi = v[1];
          if (!(v[2] >= 2. && v[2] <= 1e6) || v[2] != std::floor(v[2])) {
            throw std::runtime_error(key + ": n must be an integer >= 2! (use frequency = nu for a single frequency)");
          }
          if (!(lo > 0.) || !(hi > lo)) {
            throw std::runtime_error(key + ": need 0 < start < stop!");
          }
          const int n = (int)v[2];
          for (int i = 0; i < n; i++) {
            const double x = double(i) / (n - 1);
            out.push_back(i == n - 1 ? hi : (form == LOGSPACE ? lo * std::pow(hi / lo, x) : lo + (hi - lo) * x));
          }
        }
      }

      std::set<std::string> seen;
      for (double v : out) {
        if (!(v > 0.) || !std::isfinite(v)) {
          throw std::runtime_error("frequencies must be positive and finite");
        }
        char b[64];
        std::snprintf(b, sizeof b, "%e", v);
        if (!seen.insert(b).second) {
          throw std::runtime_error(std::string("duplicate frequency ") + b + " (they are used to name the HDF5 datasets)");
        }
      }
    }

    /**
     * @brief Parses a TOML file and places the key value pairs into a std::map<std::string, std::string>
     *
     * @param text The contents of the TOML file
     * @param kv The key value pairs of the TOML file
     * @throw std::runtime_error
     */
    inline void parse_toml(
      const std::string& text,
      std::map<std::string, std::string>& kv
    ) {
      auto trim = [](std::string t) {
        const size_t p = t.find_first_not_of(" \t\r"), q = t.find_last_not_of(" \t\r");
        return p == std::string::npos ? std::string() : t.substr(p, q - p + 1);
      };
      std::istringstream in(text);
      std::string line, section;
      int ln = 0;
      while (std::getline(in, line)) {
        ln++;
        bool q = false;  // strip a '#' comment that is not inside a string
        for (size_t i = 0; i < line.size(); i++) {
          if (line[i] == '"') {
            q = !q;
          }
          if (line[i] == '#' && !q) {
            line = line.substr(0, i);
            break;
          }
        }

        line = trim(line);
        if (line.empty()) continue;

        if (line[0] == '[') {
          if (line.back() != ']') {
            throw std::runtime_error("Line " + std::to_string(ln) + ": bad section header");
          }
          section = trim(line.substr(1, line.size() - 2));
          continue;
        }

        const size_t eq = line.find('=');
        if (eq == std::string::npos) {
          throw std::runtime_error("Line " + std::to_string(ln) + ": expected key = value");
        }

        std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
        if (val.size() >= 2 && val.front() == '"' && val.back() == '"') {
          val = val.substr(1, val.size() - 2);
        }

        std::string num;  // numbers may use TOML's digit separators
        for (char ch : val) {
          num += ch;
        }

        const std::string full = section.empty() ? key : section + "." + key;
        if (kv.count(full)) {
          throw std::runtime_error("Line " + std::to_string(ln) + ": duplicate key " + full);
        }
        kv[full] = val;
      }
    }

    /**
     * @brief Loads the config.toml file located in `filepath`, applies overrides
     *        to each one if there is a CLI argument that matches it.
     *
     * @param filepath   The config.toml file (located in input/config.toml)
     * @param overrides  The list of CLI args that override the config
     * @param config     The config struct that is written to
     * @throw std::runtime_error
     */
    inline void load_config(
      const std::string& filepath,
      const std::vector<std::string>& overrides,
      Params& params,
      core::bhac::UserParams* model_uparams
    ) {
      std::ifstream file(filepath);
      if (!file) {
        throw std::runtime_error("Cannot open file " + filepath + ". Perhaps it's a malformed TOML key-value pair?");
      }
      std::stringstream buf;
      buf << file.rdbuf();
      std::map<std::string, std::string> arg_map;

      parse_toml(buf.str(), arg_map);

      bool freq_override = false;
      for (const std::string& a : overrides) {
        const size_t eq = a.find('=');
        if (eq == std::string::npos) {
          throw std::runtime_error ("Override must be section.key=value: " + a);
        }

        // replace the whole frequency
        if (a.compare(0, 15, "camera.frequenc") == 0 && !freq_override) {
          freq_override = true;
          for (auto it = arg_map.begin(); it != arg_map.end();) {
            it = (it->first.compare(0, 15, "camera.frequenc") == 0) ?
                    arg_map.erase(it) :
                    std::next(it);
          }
        }
        arg_map[a.substr(0, eq)] = a.substr(eq + 1);
      }

      std::set<std::string> used;
      std::string missing;
      auto has = [&](const char* k) { return arg_map.count(k) > 0; };

      auto str = [&](const char* k, std::string& dst, bool required) {
        if (has(k)) {
          dst = arg_map[k];
          used.insert(k);
        } else if (required) {
          missing += std::string(" ") + k;
        }
      };

      auto dbl = [&](const char* k, double& dst, bool required) {
        std::string t;
        str(k, t, required);
        if (!t.empty()) dst = std::atof(t.c_str());
      };

      auto sclr = [&](const char* k, scalar& dst, bool required) {
        std::string t;
        str(k, t, required);
        if (!t.empty()) dst = std::atof(t.c_str());
      };

      auto integer = [&](const char* k, int& dst, bool required) {
        double v = dst;
        dbl(k, v, required);
        dst = (int)v;
      };

      auto boolean = [&](const char* k, bool& dst, bool required) {
        std::string t;
        str(k, t, required);
        if (!t.empty()) dst = (t == "true" || t == "1");
      };
      // [visual]
      boolean("visual.show_header", params.show_header, true);

      // [black_hole]
      dbl("black_hole.mass", model_uparams->mass, true);
      dbl("black_hole.distance", model_uparams->distance, true);

      // [integration]
      sclr("integration.stepsize", model_uparams->stepsize, true);
      integer("integration.max_steps", model_uparams->max_steps, true);
      dbl("integration.horizon_margin", model_uparams->horizon_margin, true);
      integer("integration.max_order", model_uparams->max_order, true);

      // [camera]
      dbl("camera.inclination", model_uparams->inclination, true);
      integer("camera.img_width", model_uparams->img_width, true);
      integer("camera.img_height", model_uparams->img_height, true);
      dbl("camera.fov_x", model_uparams->fov_x, true);
      dbl("camera.fov_y", model_uparams->fov_y, true);
      dbl("camera.r_cam", model_uparams->r_cam, true);
      dbl("camera.t_init", model_uparams->t_init, true);
      model_uparams->pa = 0.;
      model_uparams->phi_cam = 0.;
      dbl("camera.pa", model_uparams->pa, false);        // optional, deg
      dbl("camera.phi", model_uparams->phi_cam, false);  // optional, deg (needs PHX_CAMERA_PHI)
#if !PHX_CAMERA_PHI
      if (model_uparams->phi_cam != 0.) throw std::runtime_error("camera.phi needs a build with PHX_CAMERA_PHI=ON");
#endif

      try {
        std::vector<double> freqs;
        read_frequencies(arg_map, used, freqs);
        if (freqs.size() > (size_t)PHX_MAX_FREQS) {
          throw std::runtime_error("at most " + std::to_string(PHX_MAX_FREQS) + " frequencies (PHX_MAX_FREQS), got " + std::to_string(freqs.size()));
        }
        std::copy(freqs.begin(), freqs.end(), model_uparams->frequencies);
        model_uparams->nfreqs = (int)freqs.size();
      } catch (const std::exception& e) {
        throw std::runtime_error(filepath + ": " + e.what());
      }

      // [output]
      str("output.dir", params.outdir, true);
      str("output.image", params.image, false);
      integer("output.block_size", params.block_size, false);
      if (params.image != "hdf5" && params.image != "text" && params.image != "both" && params.image != "none") {
        throw std::runtime_error(filepath + ": output.image must be hdf5, text, both or none");
      }
      if (params.block_size < 1) throw std::runtime_error(filepath + ": output.block_size must be >= 1");

      // [grmhd]
      str("grmhd.dump_file", params.dump_file, true);
      str("grmhd.grid_file", params.grid_file, true);
      boolean("grmhd.gpu_decode", params.gpu_decode, true);
      integer("grmhd.nspin", model_uparams->nspin, true);
      boolean("grmhd.sfc", model_uparams->sfc, true);
      integer("grmhd.nlevels", model_uparams->nlevels, true);

      integer("grmhd.bpol", model_uparams->bpol, true);

      // [plasma]
      dbl("plasma.m_unit", model_uparams->m_unit, true);
      dbl("plasma.r_low", model_uparams->r_low, true);
      dbl("plasma.r_high", model_uparams->r_high, true);
      sclr("plasma.sigma_cut", model_uparams->sigma_cut, true);
      sclr("plasma.thetae_min", model_uparams->thetae_min, true);
      sclr("plasma.thetae_max", model_uparams->thetae_max, true);
      sclr("plasma.rt_outer", model_uparams->rt_outer, true);
      model_uparams->beta_c = 1.;
      dbl("plasma.beta_c", model_uparams->beta_c, false);  // optional, RAPTOR: 1
      if (!(model_uparams->beta_c > 0.)) throw std::runtime_error(filepath + ": plasma.beta_c must be > 0");

      if (!missing.empty()) {
        throw std::runtime_error(filepath + ": missing required key(s):" + missing);
      }

      for (const auto& e : arg_map) {
        if (!used.count(e.first)) {
          throw std::runtime_error(filepath + ": unknown key " + e.first);
        }
      }
    }
  } // namespace config
} // namespace host
