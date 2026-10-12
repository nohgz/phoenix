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

#include <string>
#include <vector>
#include <map>
#include <set>

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

      auto num = [&](const char* k, double& dst, bool required) {
        std::string t;
        str(k, t, required);
        if (!t.empty()) dst = std::atof(t.c_str());
      };

      auto integer = [&](const char* k, int& dst, bool required) {
        double v = dst;
        num(k, v, required);
        dst = (int)v;
      };

      auto boolean = [&](const char* k, bool& dst, bool required) {
        std::string t;
        str(k, t, required);
        if (!t.empty()) dst = (t == "true" || t == "1");
      };
      // [visual]
      boolean("visual.show_header", params.show_header, true);

      // [camera]
      num("camera.inclination", model_uparams->fov_x, true);
      integer("camera.img_width", model_uparams->img_width, true);
      integer("camera.img_height", model_uparams->img_height, true);
      num("camera.fov_x", model_uparams->fov_x, true);
      num("camera.fov_y", model_uparams->fov_y, true);
      num("camera.r_cam", model_uparams->r_cam, true);
      num("camera.t_init", model_uparams->t_init, true);

      try {
        read_frequencies(arg_map, used, model_uparams->frequencies);
      } catch (const std::exception& e) {
        throw std::runtime_error(filepath + ": " + e.what());
      }
      model_uparams->nfreqs = model_uparams->frequencies.size();

      // [output]
      str("output.dir", params.outdir, true);

      // [grmhd]
      str("grmhd.dump_file", params.dump_file, true);
      str("grmhd.grid_file", params.grid_file, true);
      boolean("grmhd.gpu_decode", params.gpu_decode, true);
      integer("grmhd.nspin", model_uparams->nspin, true);
      boolean("grmhd.sfc", model_uparams->sfc, true);
      integer("grmhd.nlevels", model_uparams->nlevels, true);

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
