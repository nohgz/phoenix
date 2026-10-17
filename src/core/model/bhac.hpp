#pragma once

#include <vector>

#include "core/numerics/headers.hpp"
#include "core/numerics/precision.hpp"
#include "core/geometry/metric.hpp"

namespace core {
  namespace bhac {

    // Following the RAPTOR convention for the primitive, conserved vars
    enum {
      KRHO = 0,
      UU = 1,
      U1 = 2,
      U2 = 3,
      U3 = 4,
      B1 = 5,
      B2 = 6,
      B3 = 7,
      NPRIM = 8
    };

    // Again, following RAPTORs conserved-variable slots of the dump
    enum {
      C_D = 0,
      C_S1 = 1,
      C_S2 = 2,
      C_S3 = 3,
      C_TAU = 4,
      C_B1 = 5,
      C_B2 = 6,
      C_B3 = 7,
      C_DS = 8
    };

    struct PlasmaUnits {
      double L;
      double T;
      double rho;
      double U;
      double B;
      double Ne;
    };

    /**
     * @brief Holds the grid parameters that the dump does not carry.
     */
    struct GridParams {
      int nxlone[3];
      double xprobmin[3];
      double xprobmax[3];
      double hslope;
    };

    /**
     * @brief Holds all parameters that the user can tweak.
     *        (the PHOENIX equivalent of model.in).
     *        This struct isn't v hot so I don't mind packing it poorly rn.
     */
    struct UserParams {
      // [camera]
      // The camera's inclination wrt the BH
      double inclination;
      // The image width
      int img_width;
      // The image height
      int img_height;
      // The camera horizontal field of view
      double fov_x;
      // The camera vertical field of view
      double fov_y;
      // The camera's distance to the BH
      double r_cam;
      // The time the camera sits at
      double t_init;
      // Position angle of the image on the sky (deg)
      double pa;
      // Azimuth of the camera around the BH (deg, KS phi). Only used if PHX_CAMERA_PHI (this is dumb. gonna patch it out soon.)
      double phi_cam;
      // The array of frequencies to integrate over
      double frequencies[PHX_MAX_FREQS];
      // How many frequencies are we integrating over? (indirectly set by user)
      int nfreqs;

      // [black_hole]
      // The black hole's mass (Msun)
      double mass;
      // The black hole's distance (kpc)
      double distance;

      // [integration]
      // The scaling factor for the integration step size
      scalar stepsize;
      // Maximum number of steps that each ray can take
      int max_steps;
      // used to calculate cut_in
      double horizon_margin;
      // Maximum number of times theta can wrap
      int max_order;

      // [plasma]
      // The mass unit of the simulation (g)
      double m_unit;
      // electron temperature: Rlow prescription
      double r_low;
      // electron temperature: Rhigh prescription
      double r_high;
      // plasma beta where the R_low / R_high prescription switches (RAPTOR: 1)
      double beta_c;
      // no emission where b^2/rho > sigma_cut
      scalar sigma_cut;
      // electron temperature window [m_e c^2]
      scalar thetae_min;
      scalar thetae_max;
      // no emission beyond this radius (Rg)
      scalar rt_outer;

      // [grmhd]
      // The spin's index in the dump neqpar
      int nspin;
      // The sign of the poloidal field
      int bpol;
      // The number of refinement levels
      int nlevels;
      // Should we use a Morton-ordered z curve?
      bool sfc;
    };

    struct Model {
      GridParams gparams;
      UserParams uparams;
      PlasmaUnits units;

      // from the dump header
      int nleafs, ndim, nwini, nws, cells;
      int nx[3];     // cells per block
      int ng[3];     // blocks per dimension on the base level

      // idx: [3*leaf + dim]
      scalar *lb;
      scalar *dxc;
      // idx: [(leaf*cells + cell)*NPRIM + var]
      scalar *prim;

      // for the metric
      scalar a, hslope;

      // for the plasma
      scalar gam, Ne, Bu, rhigh, rlow, ibc;  // ibc = 1 / beta_c

      // for the leaf finding table
      int level_dims[16][3];
      // leaf index per (level, i, j, k)
      // sentinel value of -1 (the thing that says go another amr level)
      int* level_map[16];
      // cell size of each level's base blocks
      scalar lvl_w[16][3];
      // simulation domain
      scalar x_lo[4], x_hi[4];

      // cosine and sine of the position angle
      double cos_pa, sin_pa;

      // for integration stuffs
      scalar r_out, cut_in;

      // for the transfer, per observing frequency
      scalar freq[PHX_MAX_FREQS];    // [Hz]
      scalar hconst[PHX_MAX_FREQS];  // m_e c^2 / (h nu) * r_g h / (m_e c^2): segment length -> optical depth
      scalar jscale[PHX_MAX_FREQS];  // nu^3: S carries the intensity directly
    };

    /**
     * @brief sqrt(det gamma_ij) of the spatial metric at (x1, x2, x3).
     */
    GRRT_SHARED_FUNCTION double get_detgamma(
      const double a,
      const double hs,
      const double x,
      const double y,
      const double z
    ) {
      const double X_u[4] = {0., x, y, z};
      double g_dd[4][4];
      core::metric::metric_dd<double>(a, hs, X_u, g_dd);
      const double det = g_dd[1][1] * (g_dd[2][2] * g_dd[3][3] - g_dd[3][2] * g_dd[2][3]) -
                         g_dd[1][2] * (g_dd[2][1] * g_dd[3][3] - g_dd[3][1] * g_dd[2][3]) +
                         g_dd[1][3] * (g_dd[2][1] * g_dd[3][2] - g_dd[2][2] * g_dd[3][1]);
      return sqrt(det);
    }

    /**
     * @brief Volume-weighted cell centroid
     *
     * @param m     Model for spin and hslope
     * @param X     Geometric cell center (x1, x2, x3)
     * @param dxc   Cell size
     * @param Xbar  Output,centroid
     */
    GRRT_SHARED_FUNCTION void coord_bar(
      const Model* model,
      const double X[3],
      const double dxc[3],
      double Xbar[3]
    ) {
      const double a = double(model->a);
      const double hslope = double(model->hslope);
      const double c1[3] = {1., 4., 1.};
      const double s[3] = {-1., 0., 1.};
      double norm = 0., xb = 0., yb = 0., zb = 0.;

      for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            const double w = c1[i] * c1[j] * c1[k];
            const double x = X[0] + dxc[0] / 2. * s[i];
            const double y = X[1] + dxc[1] / 2. * s[j];
            const double z = X[2] + dxc[2] / 2. * s[k];
            const double dg = get_detgamma(a, hslope, x, y, z);
            norm += dg * w;
            xb += dg * w * x;
            yb += dg * w * y;
            zb += dg * w * z;
          }
        }
      }

      const double vol = (dxc[0] / 6.0) * (dxc[1] / 6.0) * (dxc[2] / 6.0);
      norm *= vol;
      Xbar[0] = xb * (vol / norm);
      Xbar[1] = yb * (vol / norm);
      Xbar[2] = zb * (vol / norm);
    }

    /**
     * @brief Convert conserved variables to primitives.
     *
     * @param m     Model
     * @param lb3   Block origin of the leaf
     * @param dxc3  Cell size of the leaf
     * @param cons  The leaf's raw block, nwini arrays of `cells` doubles
     * @param cell  Cell index within the leaf
     * @param out   Output, NPRIM values (converted to scalar on the way out)
     */
    GRRT_SHARED_FUNCTION void convert2prim(
      const Model* model,
      const scalar lb3[3],
      const scalar dxc3[3],
      const double* cons,
      int cell,
      scalar* prim  // NPRIM outputs
    ) {

// Cv = Conserved variable. This is just for convenience only here
#define Cv(nw) (cons[(size_t)(nw)*cells + cell])

      const int cells = model->cells;
      const int LFAC = model->nwini - 2, XI = model->nwini - 1;
      const int loc[3] = {
        cell % model->nx[0],
       (cell / model->nx[0]) % model->nx[1],
        cell / (model->nx[0] * model->nx[1])
      };

      double X[3], dxc[3];
      for (int dim = 0; dim < 3; dim++) {
        dxc[dim] = (double)dxc3[dim];
        X[dim] = (double)lb3[dim] + (loc[dim] + 0.5) * (double)dxc3[dim];
      }

      for (int var = 0; var < NPRIM; var++) {
        prim[var] = scalar(0);
      }

      if (exp(X[0]) > 1.0) {
        double Xbar[3];
        coord_bar(model, X, dxc, Xbar);
        const double X_u[4] = {0., Xbar[0], Xbar[1], Xbar[2]};

        double g_dd[4][4], g_uu[4][4];
        core::metric::metric_dd_uu<double>(
          double(model->a),
          double(model->hslope),
          X_u,
          g_dd,
          g_uu
        );

        const double BS = Cv(C_S1) * Cv(C_B1)
                        + Cv(C_S2) * Cv(C_B2)
                        + Cv(C_S3) * Cv(C_B3);

        double Bsq = 0.;
        double B_d[4] = {0, 0, 0, 0};
        double S_u[4] = {0, 0, 0, 0};
        double gamma_uu[4][4];

        for (int i = 1; i < 4; i++) {
          for (int j = 1; j < 4; j++) {
            gamma_uu[i][j] = g_uu[i][j] + g_uu[0][i] * g_uu[0][j] / (-g_uu[0][0]);
          }
        }

        for (int j = 1; j < 4; j++) {
          for (int i = 1; i < 4; i++) {
            S_u[j] += gamma_uu[i][j] * Cv(C_S1 + i - 1);
            B_d[j] += g_dd[i][j] * Cv(C_B1 + i - 1);
          }
        }

        for (int i = 1; i < 4; i++) {
          Bsq += B_d[i] * Cv(C_B1 + i - 1);
        }

        const double gam = double(model->gam);
        const double lfac = Cv(LFAC), xi = Cv(XI);

        const double rho = Cv(C_D) / lfac;
        double uu = (gam - 1.) / gam * (xi / pow(lfac, 2.) - rho) / (gam - 1.);

        if (uu < 0) {  // entropy fallback
          uu = (Cv(C_DS) / Cv(C_D)) * pow(rho, gam - 1) / (gam - 1.);
          if (uu < 0) uu = 1e-14;
        }

        const double den = xi + Bsq;
        const double u1 = (S_u[1] / den + Cv(C_B1) * BS / (xi * den)) * lfac;
        const double u2 = (S_u[2] / den + Cv(C_B2) * BS / (xi * den)) * lfac;
        const double u3 = (S_u[3] / den + Cv(C_B3) * BS / (xi * den)) * lfac;
        const double b1 = model->uparams.bpol * Cv(C_B1);
        const double b2 = model->uparams.bpol * Cv(C_B2);
        const double b3 = model->uparams.bpol * Cv(C_B3);

        prim[KRHO] = scalar(rho);
        prim[UU]   = scalar(uu);
        prim[U1]   = scalar(u1);
        prim[U2]   = scalar(u2);
        prim[U3]   = scalar(u3);
        prim[B1]   = scalar(b1);
        prim[B2]   = scalar(b2);
        prim[B3]   = scalar(b3);
      }
#undef Cv
    }
  }
}