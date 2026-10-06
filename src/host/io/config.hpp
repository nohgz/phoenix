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

namespace host {
  namespace config {
    /**
    * @brief See the comments in the config.toml file for descriptions of these
    *        variables.
    */
    struct Config {
      // [visual]
      bool show_header;
    };

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
      Config& config
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

      boolean("visual.show_header", config.show_header, true);

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
