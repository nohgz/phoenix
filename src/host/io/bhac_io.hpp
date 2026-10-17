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
#include "host/formatting/fancy_printing.hpp"
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
      core::bhac::GridParams* gparams,
      const std::string& gridfile
    ) {
      std::ifstream grid_stream(gridfile);
      if (!grid_stream) {
        throw std::runtime_error("cannot open grid file " + gridfile);
      }
      std::vector<double> gval;
      std::string line;
      while (std::getline(grid_stream, line)) {
        std::istringstream ss(line);
        std::string tok, last;
        while (ss >> tok) last = tok;
        if (!last.empty()) gval.push_back(std::atof(last.c_str()));
      }
      if (gval.size() < 10) {
        throw std::runtime_error(gridfile + ": expected 10 lines (nxlone x3, xprobmin x3, xprobmax x3, hslope)");
      }
      for (int d = 0; d < 3; d++) {
        gparams->nxlone[d] = (int)gval[d];
        gparams->xprobmin[d] = gval[3 + d];
        gparams->xprobmax[d] = gval[6 + d];
      }

      gparams->hslope = gval[9];
      if (gparams->nxlone[0] <= 0 || gparams->nxlone[1] <= 0 || gparams->nxlone[2] <= 0 ||
          gparams->xprobmax[0] <= gparams->xprobmin[0] || gparams->xprobmax[1] <= gparams->xprobmin[1] ||
          gparams->xprobmax[2] <= gparams->xprobmin[2]) {
        throw std::runtime_error(gridfile + ": nxlone / xprobmin / xprobmax not set (or empty domain)");
      }
    }

    /**
     * @brief Finalize the model with info that the dump's tables dont have
     *
     * @param model   The model
     * @param neqpar  The dump's equation parameters
     * @throw std::runtime_error
     */
    static void finalize_model(
      core::bhac::Model* model,
      const std::vector<double>& neqpar
    ) {
      namespace k = core::constants;
      const core::bhac::UserParams& uparams = model->uparams;

      // -- metric and equation of state, from the dump
      const double spin = neqpar[uparams.nspin];
      if (std::fabs(spin) > 1.) {
        throw std::runtime_error("Spin " + std::to_string(spin) + " unphysical (wrong grmhd.nspin?)");
      }
      model->a = scalar(spin);
      model->hslope = scalar(model->gparams.hslope);
      model->gam = scalar(neqpar[0]);

      // mirrors RAPTOR's BHAC model.c:649-655
      for (int dim = 1; dim < 3; dim++) {
        model->gparams.xprobmin[dim] *= 2. * M_PI;
        model->gparams.xprobmax[dim] *= 2. * M_PI;
      }
      model->cells = model->nx[0] * model->nx[1] * model->nx[2];
      for (int dim = 0; dim < 3; dim++) {
        if (model->gparams.nxlone[dim] % model->nx[dim] != 0) {
          throw std::runtime_error("nxlone[" + std::to_string(dim) + "] = " + std::to_string(model->gparams.nxlone[dim]) +
                                    " is not a multiple of the dump's block size " + std::to_string(model->nx[dim]) +
                                    ": wrong grid file");
        }
        model->ng[dim] = model->gparams.nxlone[dim] / model->nx[dim];
        model->x_lo[1 + dim] = scalar(model->gparams.xprobmin[dim]);
        model->x_hi[1 + dim] = scalar(model->gparams.xprobmax[dim]);
      }

      for (int lvl = 0; lvl < 16; lvl++) {
        for (int dim = 0; dim < 3; dim++) {
          model->lvl_w[lvl][dim] = scalar(
            (model->gparams.xprobmax[dim] - model->gparams.xprobmin[dim]) / model->ng[dim] / (double)(1 << lvl)
          );
        }
      }

      // units and plasma scaling
      model->units.L = k::G_GRAV * uparams.mass * k::M_SUN / k::SPEED_OF_LIGHT / k::SPEED_OF_LIGHT;
      model->units.T = model->units.L / k::SPEED_OF_LIGHT;
      const double rho = model->uparams.m_unit / std::pow(model->units.L, 3);
      model->units.rho = rho;
      model->units.B = core::constants::SPEED_OF_LIGHT * std::sqrt(4. * M_PI * rho);
      model->units.Ne = rho / (core::constants::M_PROTON + core::constants::M_ELECTRON);
      model->Bu = scalar(model->units.B);
      model->Ne = scalar(model->units.Ne);
      model->rhigh = scalar(model->uparams.r_high);
      model->rlow = scalar(model->uparams.r_low);
      model->ibc = scalar(1. / model->uparams.beta_c);
      model->units.U = rho * k::SPEED_OF_LIGHT * k::SPEED_OF_LIGHT;

      // camera and integration
      model->cos_pa = std::cos(uparams.pa / 180. * M_PI);
      model->sin_pa = std::sin(uparams.pa / 180. * M_PI);
      model->cut_in = scalar((1. + std::sqrt(1. - spin * spin)) * (1. + uparams.horizon_margin));
      model->r_out = scalar(1.1 * uparams.r_cam);

      const double c_const = model->units.L * k::H_PLANCK / (k::M_ELECTRON * k::SPEED_OF_LIGHT * k::SPEED_OF_LIGHT);
      for (int f = 0; f < uparams.nfreqs; f++) {
        const double nu = uparams.frequencies[f];
        model->freq[f] = scalar(nu);
        model->hconst[f] = scalar((k::M_ELECTRON * k::SPEED_OF_LIGHT * k::SPEED_OF_LIGHT) / (k::H_PLANCK * nu) * c_const);
        model->jscale[f] = scalar(nu * nu * nu);
      }
    }

    /**
     * @brief Leaf index per triplet (level, block index).
     *        Yields a O(1) block lookup for the sampler.
     */
    static void build_level_map(
      core::bhac::Model* model,
      const std::vector<int>& level,
      const std::vector<int>& ind
    ) {
      if (model->uparams.nlevels == -1) {
        // O(n) to get the maximum refinement level if not provided by the user
        for (int leaf = 0; leaf < model->nleafs; leaf++) {
          model->uparams.nlevels = std::max(model->uparams.nlevels, level[leaf]);
        }
        if (model->uparams.nlevels > 16) {
          throw std::runtime_error(std::to_string(model->uparams.nlevels) + " AMR levels unsupported (max 16)");
        }
      }

      for (int lvl = 0; lvl < model->uparams.nlevels; lvl++) {
        const int sc = 1 << lvl;
        for (int dim = 0; dim < 3; dim++) {
          model->level_dims[lvl][dim] = model->ng[dim] * sc;
        }

        const int total =
            model->level_dims[lvl][0]
          * model->level_dims[lvl][1]
          * model->level_dims[lvl][2];
        model->level_map[lvl] = (int*)std::malloc(total * sizeof(int));
        for (int i = 0; i < total; i++) {
          model->level_map[lvl][i] = -1;
        }

        for (int leaf = 0; leaf < model->nleafs; leaf++) {
          if (level[leaf] - 1 != lvl) continue;

          const int i = ind[3 * leaf + 0];
          const int j = ind[3 * leaf + 1];
          const int k = ind[3 * leaf + 2];

          if (
            i < 0 || i >= model->level_dims[lvl][0] ||
            j < 0 || j >= model->level_dims[lvl][1] ||
            k < 0 || k >= model->level_dims[lvl][2]
          ) {
            host::fancy_printing::infomsg("Leaf %d index (%d,%d,%d) outside level-%d dims\n", leaf, i, j, k, lvl + 1);
            continue;
          }
          const int idx = (i * model->level_dims[lvl][1] + j) * model->level_dims[lvl][2] + k;
          if (model->level_map[lvl][idx] != -1) {
            host::fancy_printing::failmsg("Level-map collision at leaf %d\n", leaf);
          }
          model->level_map[lvl][idx] = leaf;
        }
      }
    }

    /**
     * @brief Frees the host tables of a Model.
     */
    static void free_model(core::bhac::Model* model) {
      std::free(model->lb);
      std::free(model->dxc);
      std::free(model->prim);
      for (int lvl = 0; lvl < 16; lvl++) {
        std::free(model->level_map[lvl]);
      }
      std::memset(model, 0, sizeof *model);
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
      const std::string& dump_filepath,
      bool gpu_decode
    ) {
      FILE *dumpfile = std::fopen(dump_filepath.c_str(), "rb");

      if (!dumpfile) {
        throw std::runtime_error("Can't open " + dump_filepath);
      }

      try {
        // unused integer values are dumped here
        int idump;
        // unused double values are dumped here
        double ddump;

        int neqparini;

        // -- header lives at the end of the file
        std::fseek(dumpfile, 0, SEEK_END);
        long offset = -40;
        std::fseek(dumpfile, offset, SEEK_CUR);

        bool read_ok = read_int(dumpfile, &model->nleafs) && read_int(dumpfile, &idump)
          && read_int(dumpfile, &model->ndim) && read_int(dumpfile, &idump)
          && read_int(dumpfile, &model->nwini) && read_int(dumpfile, &model->nws)
          && read_int(dumpfile, &neqparini) && read_int(dumpfile, &idump)
          && read_dbl(dumpfile, &ddump);

        if (!read_ok) {
          throw std::runtime_error("Bad read in the dump header. Is " + dump_filepath + " a BHAC dump?");
        }
        if (model->ndim != 3) {
          throw std::runtime_error("Only 3D dumps are supported! (ndim = " + std::to_string(model->ndim) + ")");
        }

        offset = offset - (model->ndim * 4 + neqparini * 8);
        std::fseek(dumpfile, offset, SEEK_CUR);
        for (int dim = 0; dim < model->ndim; dim++) {
          read_ok = read_ok && read_int(dumpfile, &model->nx[dim]);
        }

        std::vector<double> neqpar(neqparini);
        for (int i = 0; i < neqparini; i++) {
          read_ok = read_ok && read_dbl(dumpfile, &neqpar[i]);
        }

        if (!read_ok || model->uparams.nspin < 0 || model->uparams.nspin >= neqparini) {
          throw std::runtime_error("Short read in the dump header! (wrong grmhd.nspin?)");
        }

        // everything derived from the parameters and the header, before anything is decoded
        finalize_model(model, neqpar);

        // -- loraxmaxxing (walking the tree for leaves)
        const long size_blk = (long)model->cells * model->nwini * 8;
        const long stag_per_leaf = (long)(model->nx[0] + 1) * (model->nx[1] + 1) * (model->nx[2] + 1) * model->nws * 8;
        std::fseek(dumpfile, (long)model->nleafs * (size_blk + stag_per_leaf), SEEK_SET);
        std::vector<int> ind, level;

        if (model->uparams.sfc) {
          std::vector<std::pair<uint64_t, int>> blocks;  // (Morton code, i*ng1*ng2 + j*ng2 + k)
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
              dumpfile,
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
                read_node(dumpfile, model->ndim, 1, i, j, k, ind, level, read_ok);
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

        // -- block geometry
        const size_t n3 = 3 * (size_t)model->nleafs;
        model->lb = (scalar*)std::malloc(n3 * sizeof(scalar));
        model->dxc = (scalar*)std::malloc(n3 * sizeof(scalar));
        double dx1[3], dxc0[3];

        for (int dim = 0; dim < 3; dim++) {
          dx1[dim] =  (model->gparams.xprobmax[dim] - model->gparams.xprobmin[dim]) / model->ng[dim];
          dxc0[dim] = (model->gparams.xprobmax[dim] - model->gparams.xprobmin[dim]) / model->gparams.nxlone[dim];
        }

        for (int leaf = 0; leaf < model->nleafs; leaf++) {
          const double inv_scale = 1.0 / (double)(1 << (level[leaf] - 1));
          for (int dim = 0; dim < 3; dim++) {
            model->lb[3 * leaf + dim] = scalar(model->gparams.xprobmin[dim] + ind[3 * leaf + dim] * dx1[dim] * inv_scale);
            model->dxc[3 * leaf + dim] = scalar(dxc0[dim] * inv_scale);
          }
        }
        build_level_map(model, level, ind);

        if (!gpu_decode) {
          // Two phases per chunk of leaves:
          // 1. serial file reads into a buffer
          // 2. conserved -> primitive over the chunk in parallel
          model->prim = (scalar*)std::calloc((size_t)model->nleafs * model->cells * core::bhac::NPRIM, sizeof(scalar));
          const int CHUNK = 256;
          const size_t leaf_doubles = (size_t)model->nwini * model->cells;
          std::vector<double> raw((size_t)CHUNK * leaf_doubles);
          std::fseek(dumpfile, 0, SEEK_SET);
          for (int leaf0 = 0; leaf0 < model->nleafs; leaf0 += CHUNK) {
            const int nl = std::min(CHUNK, model->nleafs - leaf0);
            // **** Phase 1 ****
            for (int q = 0; q < nl; q++) {
              if (std::fread(raw.data() + (size_t)q * leaf_doubles, sizeof(double), leaf_doubles, dumpfile) != leaf_doubles) {
                throw std::runtime_error("short read in the data, leaf " + std::to_string(leaf0 + q));
              }
              std::fseek(dumpfile, stag_per_leaf, SEEK_CUR);
            }

            // **** Phase 2 ****
#pragma omp parallel for schedule(dynamic, 1)
            for (int q = 0; q < nl; q++) {
              const int leaf = leaf0 + q;
              const double* leaf_buf = raw.data() + (size_t)q * leaf_doubles;
              for (int cell = 0; cell < model->cells; cell++) {
                core::bhac::convert2prim(
                  model,
                  &model->lb[3 * leaf],
                  &model->dxc[3 * leaf],
                  leaf_buf,
                  cell,
                  &model->prim[((size_t)leaf * model->cells + cell) * core::bhac::NPRIM]
                );
              }
            }
          }
        }
      } catch (...) {
        std::fclose(dumpfile);
        free_model(model);
        throw;
      }
      std::fclose(dumpfile);
    }
  } // namespace bhac_io
} // namespace host
