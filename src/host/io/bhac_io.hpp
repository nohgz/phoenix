/**
 * @file bhac_io.hpp
 * @author C. Noga & J. Davelaar (github.com/nohgz)
 * @brief Contains the functions necessary to load BHAC data dumps
 * @version 1.0
 * @date 2026-10-06
 * 
 * GNU GPL v3.0
 * 
 */

#pragma once

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <stdexcept>

#include "core/model/headers.hpp"
#include "host/io/config.hpp"

namespace host {
  namespace bhac_io {
    /**
     * @brief 
     * @author J. Davelaar
     * 
     * @param dims 
     * @param i 
     * @param new_i 
     * @param new_j 
     * @param new_k 
     * @param ip 
     * @param jp 
     * @param kp 
     */
    static void new_index(
      int dims,
      int i,
      int *new_i,
      int *new_j,
      int *new_k,
      int ip,
      int jp,
      int kp
    ) {
      if (dims == 1) {
        *new_i = 2 * (ip) + i;
        *new_j = 0;
        *new_k = 0;
      }

      if (dims == 2) {
        *new_i = 2 * (ip) + i % 2;
        *new_j = 2 * (jp) + i / 2;
        *new_k = 0;
      }

      if (dims == 3) {
        *new_i = 2 * (ip) + i % 2;
        *new_j = 2 * (jp) + ((int)(i / 2.)) % 2;
        *new_k = 2 * (kp) + (int)(i / 4.);
      }
    }

    /**
     * @brief Reads one node of the AMR tree. Leaves are appended to ind 
     *        (3 per leaf) and lev, in the order their data appears in the file.
     * @author J. Davelaar, adapted to C++ by C. Noga
     * 
     * @param f 
     * @param ndim 
     * @param level 
     * @param i 
     * @param j 
     * @param k 
     * @param ind 
     * @param lev 
     * @param ok 
     */
    static void read_node(
      FILE* f,
      int ndim,
      int level,
      int i,
      int j,
      int k,
      std::vector<int>& ind,
      std::vector<int>& lev,
      bool& ok
    ) {
      if (!ok) return;
      int is_leaf = 0;
      if (std::fread(&is_leaf, sizeof(int), 1, f) != 1) {
        ok = false;
        return;
      }
      if (is_leaf) {
        ind.push_back(i);
        ind.push_back(j);
        ind.push_back(k);
        lev.push_back(level);
        return;
      }
      for (int ich = 0; ich < (1 << ndim); ich++) {
        int ci = 0, cj = 0, ck = 0;
        new_index(ndim, ich, &ci, &cj, &ck, i, j, k);
        read_node(f, ndim, level + 1, ci, cj, ck, ind, lev, ok);
      }
    }

    /**
     * @brief Returns the morton code which interleaves bits of a, b, and c
     * @author J. Davelaar
     * 
     * @param a 
     * @param b 
     * @param c 
     * @return uint64_t 
     */
    static inline uint64_t morton3(unsigned a, unsigned b, unsigned c) {
      uint64_t ans = 0;
      for (uint64_t i = 0; i < (sizeof(uint64_t) * 8) / 3; ++i) {
        ans |= (((uint64_t)a & ((uint64_t)1 << i)) << (2 * i)) | (((uint64_t)b & ((uint64_t)1 << i)) << (2 * i + 1)) |
               (((uint64_t)c & ((uint64_t)1 << i)) << (2 * i + 2));
      }
      return ans;
    }

    static inline bool read_int(FILE* f, int* v) {
      return std::fread(v, sizeof(int), 1, f) == 1;
    }
    static inline bool read_dbl(FILE* f, double* v) {
      return std::fread(v, sizeof(double), 1, f) == 1;
    }

    /**
     * @brief Reads a BHAC base-grid file 
     * @throw std::runtime_error
     */
    inline void load_bhac_grid(
      core::bhac::GridParams& gparams,
      std::string gridfile
    ) {
      std::ifstream gi(gridfile);
      if (!gi) {
        throw std::runtime_error("cannot open grid file " + gridfile);
      }
      std::vector<double> gv;
      std::string line;
      while (std::getline(gi, line)) {
        std::istringstream ss(line);
        std::string tok, last;
        while (ss >> tok) last = tok;
        if (!last.empty()) gv.push_back(std::atof(last.c_str()));
      }
      if (gv.size() < 10) {
        throw std::runtime_error(gridfile + ": expected 10 lines (nxlone x3, xprobmin x3, xprobmax x3, hslope)");
      }
      for (int d = 0; d < 3; d++) {
        gparams.nxlone[d] = (int)gv[d];
        gparams.xprobmin[d] = gv[3 + d];
        gparams.xprobmax[d] = gv[6 + d];
      }
      gparams.hslope = gv[9];
    }

    /**
     * @brief Prepares the Model struct for rendering
     * 
     * @param model The information contained in the run's BHAC model
     * @param params The config.toml parameters
     * @param uparams The bhac-specific 
     */
    static void load(
      core::bhac::Model* model,
      std::string dump_filepath
    ) {
      std::memset(model, 0, sizeof *model);
      FILE *infile = std::fopen(dump_filepath.c_str(), "rb");

      if (!infile) {
        throw std::runtime_error("Can't open " + dump_filepath);
      }

      try {
        // unused integer values are dumped here
        int idump;
        // unused double values are dumped here
        double ddump;

        int neqparini;

        std::fseek(infile, 0, SEEK_END);
        long offset = -40;
        std::fseek(infile, offset, SEEK_CUR);

        bool read_ok = read_int(infile, &model->nleafs) && read_int(infile, &idump) 
          && read_int(infile, &model->ndim) && read_int(infile, &idump)
          && read_int(infile, &model->nwini) && read_int(infile, &model->nws)
          && read_int(infile, &neqparini) && read_int(infile, &idump)
          && read_dbl(infile, &ddump);

        if (!read_ok) {
          throw std::runtime_error("Bad read in the dump header. Is " + dump_filepath + " a BHAC dump?");
        }
        if (model->ndim != 3) {
          throw std::runtime_error("Only 3D dumps are supported! (ndim = " + std::to_string(model->ndim) + ")");
        }

        offset = offset - (model->ndim * 4 + neqparini * 8);
        std::fseek(infile, offset, SEEK_CUR);
        for (int dim = 0; dim < model->ndim; dim++) {
          read_ok = read_ok && read_int(infile, &model->nx[dim]);
        }

        std::vector<double> neqpar(neqparini);
        for (int dim = 0; dim < model->ndim; dim++) {
          read_ok = read_ok && read_dbl(infile, &neqpar[dim]);
        }

        if (!read_ok || model->uparams.nspin >= neqparini) {
          throw std::runtime_error("Short read in the dump header! (wrong grmhd.nspin?)");
        }

        for (int dim = 0; dim < 3; dim++) {
          if (model->gparams.nxlone[dim] <= 0 || 
              model->gparams.xprobmax[dim] <= model->gparams.xprobmin[dim]) {
            throw std::runtime_error("nxlone / xprobmin / xprobmax not set (or empty domain) in the grid file!");
          }
        }

        // mirrors RAPTOR's BHAC model.c:649-655
        for (int dim = 1; dim < 3; dim++) {
          model->gparams.xprobmin[dim] *= 2. * M_PI;
          model->gparams.xprobmax[dim] *= 2. * M_PI;
        }

        model->cells = model->nx[0] * model->nx[1] * model->nx[2];

        for (int dim = 0; dim < 3; dim++) {
          if (model->nxlone[dim] % model->nx[dim] != 0) {
            throw std::runtime_error("nxlone[" + std::to_string(dim) + "] = " + std::to_string(model->nxlone[dim]) +
                                     " is not a multiple of the dump's block size " + std::to_string(model->nx[dim]) +
                                     ": wrong grid file");
          }
          model->ng[dim] = model->nxlone[dim] / model->nx[dim];
        }

        // loraxmaxxing (walking the tree for leaves)
        const long size_blk = (long)model->cells * model->nwini * 8;
        const long stag_per_leaf = (long)(model->nx[0] + 1) * (model->nx[1] + 1) * (model->nx[2] + 1) * model->nws * 8;
        std::fseek(infile, (long)model->nleafs * (size_blk + stag_per_leaf), SEEK_SET);
        std::vector<int> ind, level;
        
        if (model->uparams.sfc) {
          std::vector<std::pair<uint64_t, int> > blocks;  // (Morton code, i*ng1*ng2 + j*ng2 + k)
          for (int i = 0; i < model->ng[0]; i++) {
            for (int j = 0; j < model->ng[1]; j++) {
              for (int k = 0; k < model->ng[2]; k++) {
                blocks.push_back(std::make_pair(morton3(i, j, k), (i * model->ng[1] + j) * model->ng[2] +  k));
              }
            }
          }

          std::sort(blocks.begin(), blocks.end());

          for (size_t b = 0; b < blocks.size(); b++) {
            const int n = blocks[b].second;
            read_node(
              infile,
              model->ndim, 1,
              n / (model->ng[1] * model->ng[2]),
              (n / model->ng[2]) % model->ng[1],
              n % model->ng[2],
              ind,
              level,
              read_ok
            );
          }
        } else {
          for (int k = 0; k < model->ng[2]; k++) {
            for (int j = 0; j < model->ng[1]; j++) {
              for (int i = 0; i < model->ng[0]; i++) {
                read_node(infile, model->ndim, 1, i, j, k, ind, level, read_ok);
              }
            }
          }
        }

        if (!read_ok) {
          throw std::runtime_error("short read in the tree (wrong nxlone / grmhd.sfc for this dump?)");
        }
        
        if ((int)level.size() != model->nleafs) {
          throw std::runtime_error("grid mismatch: tree has " + std::to_string(level.size()) + " leaves, header says " +
                                   std::to_string(model->nleafs) + " (grmhd.sfc / grid file?)");
        }
      } catch (...) {

      }
    }
  }
}
