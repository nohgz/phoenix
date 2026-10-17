/**
 * @file output.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Image / spectrum / HDF5 output.
 *
 * GNU GPL v3.0
 *
 * Results are stored as res[(i * img_height + j) * nfreq + f]: pixel (i, j) = alpha index i, beta index j, one
 * FreqOut per observing frequency.
 */

#pragma once

#include <hdf5.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "core/geometry/headers.hpp"
#include "core/render/headers.hpp"

namespace host {
  namespace output {
    using core::bhac::Model;
    using core::pixel::FreqOut;
    using Results = std::vector<FreqOut>;

    static inline const FreqOut& result_at(const Model* model, const Results& res, int i, int j, int f) {
      return res[((size_t)i * model->uparams.img_height + j) * model->uparams.nfreqs + f];
    }

    static inline void impact(const Model* model, int i, int j, double* alpha, double* beta) {
      const core::bhac::UserParams& up = model->uparams;
      core::camera::pixel_impact(i, j, up.img_width, up.img_height, up.fov_x, up.fov_y, alpha, beta);
    }

    /**
     * @brief Integrated Stokes flux of frequency index `fi` [Jy]. Pixels with a non-finite Stokes value are skipped
     *        (the driver reports how many there are).
     */
    static inline void compute_spec(const Model* model, const Results& res, int fi, double S[4]) {
      const double dOmega = core::camera::pixel_solid_angle(model);
      for (int q = 0; q < 4; q++) S[q] = 0.;
      for (int i = 0; i < model->uparams.img_width; i++)
        for (int j = 0; j < model->uparams.img_height; j++) {
          const FreqOut& r = result_at(model, res, i, j, fi);
          if (!std::isfinite(r.IQUV[0] + r.IQUV[1] + r.IQUV[2] + r.IQUV[3])) continue;
          for (int q = 0; q < 4; q++) S[q] += r.IQUV[q] * dOmega;
        }
      for (int q = 0; q < 4; q++) S[q] *= core::constants::JANSKY;
    }

    /**
     * @brief Mean optical depths over the finite pixels (diagnostics): tau (absorption), tauF (|Faraday rotation|).
     */
    static inline void compute_tau(const Model* model, const Results& res, int fi, double T[2]) {
      double s0 = 0., s1 = 0.;
      size_t n = 0;
      for (int i = 0; i < model->uparams.img_width; i++)
        for (int j = 0; j < model->uparams.img_height; j++) {
          const FreqOut& r = result_at(model, res, i, j, fi);
          if (!std::isfinite(r.IQUV[0] + r.IQUV[1] + r.IQUV[2] + r.IQUV[3] + r.tau + r.tauF)) continue;
          s0 += r.tau;
          s1 += r.tauF;
          n++;
        }
      T[0] = n ? s0 / (double)n : 0.;
      T[1] = n ? s1 / (double)n : 0.;
    }

    /**
     * @brief Text image of frequency index `fi`, one row per pixel (i outer, j inner): x, y [arcsec], I, Q, U, V
     *        [Jy / pixel], tau, tauF.
     */
    static inline bool write_image(const std::string& path, const Model* model, const Results& res, int fi) {
      FILE* f = std::fopen(path.c_str(), "w");
      if (!f) return false;
      const double unit = core::constants::JANSKY * core::camera::pixel_solid_angle(model);
      const double arcsec = 206265.0 * model->units.L / (model->uparams.distance * core::constants::KPC);
      for (int i = 0; i < model->uparams.img_width; i++)
        for (int j = 0; j < model->uparams.img_height; j++) {
          double a, b;
          impact(model, i, j, &a, &b);
          const FreqOut& r = result_at(model, res, i, j, fi);
          std::fprintf(f, "%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n", a * arcsec, b * arcsec,
                       r.IQUV[0] * unit, r.IQUV[1] * unit, r.IQUV[2] * unit, r.IQUV[3] * unit, r.tau, r.tauF);
        }
      std::fclose(f);
      return true;
    }

    /**
     * @brief One row per frequency: nu I Q U V [Jy], mean tau, mean tauF.
     */
    static inline bool write_spectrum(const std::string& path, const Model* model, const Results& res) {
      FILE* f = std::fopen(path.c_str(), "w");
      if (!f) return false;
      for (size_t fi = 0; fi < (size_t)model->uparams.nfreqs; fi++) {
        double S[4], T[2];
        compute_spec(model, res, (int)fi, S);
        compute_tau(model, res, (int)fi, T);
        std::fprintf(f, "%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n", model->uparams.frequencies[fi], S[0], S[1], S[2],
                     S[3], T[0], T[1]);
      }
      std::fclose(f);
      return true;
    }

    // -- RAPTOR-format HDF5 image: file img_data_<int t_init>.h5 holding 2-D double datasets [tot_blocks][bs*bs]
    // named I%e, tau%e, Q%e, U%e, V%e, tauF%e (%e = frequency in Hz, one set per frequency) plus alpha, beta (Rg).
    // I, Q, U, V are specific intensity in Jy/sr (I_cgs * 1e23; NOT Jy per pixel, unlike image_cpu.dat).
    // Layout (level-1 blocks, no camera refinement): block = bi * nblocks + bj, pixel = pi * bs + pj, image pixel
    // (i, j) = (bi * bs + pi, bj * bs + pj); alpha follows i, beta follows j.

    enum { H5_I, H5_Q, H5_U, H5_V, H5_TAU, H5_TAUF, H5_ALPHA, H5_BETA };

    static inline double h5_value(const Model* model, const Results& res, int what, int f, int i, int j) {
      if (what <= H5_V) return result_at(model, res, i, j, f).IQUV[what] * core::constants::JANSKY;
      if (what == H5_TAU) return result_at(model, res, i, j, f).tau;
      if (what == H5_TAUF) return result_at(model, res, i, j, f).tauF;
      double a, b;
      impact(model, i, j, &a, &b);
      return what == H5_ALPHA ? a : b;
    }

    static inline void h5_put(hid_t file, const std::string& name, const Model* model, int bs, const Results& res,
                              int what, int f, std::vector<double>& buf) {
      const int nbj = model->uparams.img_height / bs, nblocks = (model->uparams.img_width / bs) * nbj, npb = bs * bs;
      for (int b = 0; b < nblocks; b++)
        for (int p = 0; p < npb; p++) {
          const int i = (b / nbj) * bs + p / bs, j = (b % nbj) * bs + p % bs;
          buf[(size_t)b * npb + p] = h5_value(model, res, what, f, i, j);
        }
      hsize_t dims[2] = {(hsize_t)nblocks, (hsize_t)npb};
      hid_t sp = H5Screate_simple(2, dims, NULL);
      hid_t ds = H5Dcreate2(file, name.c_str(), H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      H5Dclose(ds);
      H5Sclose(sp);
    }

    static inline bool write_image_hdf5(const std::string& path, const Model* model, int bs, const Results& res) {
      if (model->uparams.img_width % bs || model->uparams.img_height % bs) {
        std::fprintf(stderr, "output.block_size %d must divide the image size %dx%d\n", bs, model->uparams.img_width, model->uparams.img_height);
        return false;
      }
      hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (file < 0) return false;
      std::vector<double> buf((size_t)(model->uparams.img_width / bs) * (model->uparams.img_height / bs) * bs * bs);
      for (size_t fi = 0; fi < (size_t)model->uparams.nfreqs; fi++) {
        char fb[64];
        std::snprintf(fb, sizeof fb, "%e", model->uparams.frequencies[fi]);
        const std::string fs = fb;
        h5_put(file, "I" + fs, model, bs, res, H5_I, (int)fi, buf);
        h5_put(file, "tau" + fs, model, bs, res, H5_TAU, (int)fi, buf);
        h5_put(file, "Q" + fs, model, bs, res, H5_Q, (int)fi, buf);
        h5_put(file, "U" + fs, model, bs, res, H5_U, (int)fi, buf);
        h5_put(file, "V" + fs, model, bs, res, H5_V, (int)fi, buf);
        h5_put(file, "tauF" + fs, model, bs, res, H5_TAUF, (int)fi, buf);
      }
      h5_put(file, "alpha", model, bs, res, H5_ALPHA, 0, buf);
      h5_put(file, "beta", model, bs, res, H5_BETA, 0, buf);
      H5Fclose(file);
      return true;
    }
  }
}
