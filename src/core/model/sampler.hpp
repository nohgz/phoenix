#pragma once

#include "core/numerics/headers.hpp"
#include "core/model/bhac.hpp"

namespace core {
  namespace sampler {

    /**
     * @brief Fluid at a point for the transfer stage
     */
    struct Fluid {
      // plasma four velocity (code units)
      scalar U_u[4];
      // magnetic four vector (code units)
      scalar B_u[4];
      // electron density
      scalar n_e;
      scalar B;
      scalar theta_e;
      // the four position actually sampled
      scalar X_u[4];
    };

    /**
     * @brief Wraps the four-position in place. The reflected phi (+pi) is never re-wrapped, so
     *        phi can reach 3 pi and then misses the grid.
     */
    GRRT_SHARED_FUNCTION void wrap_coords(scalar x[4]) {
      const scalar TP = scalar(2. * M_PI), PI = scalar(M_PI);
      x[3] = fmod_t(x[3], TP);
      x[2] = fmod_t(x[2], PI) - scalar(1e-6);
      if (x[3] < scalar(0)) x[3] = TP + x[3];
      if (x[2] < scalar(0)) {
        x[2] = -x[2];
        x[3] = PI + x[3];
      }
    }

    /**
     * @brief Is this position currently outside of the simulation domain?
     */
    GRRT_SHARED_FUNCTION int outside_domain(
      const bhac::Model* model,
      const scalar x[4]
    ) {
      return x[1] > model->x_hi[1] || x[1] < model->x_lo[1] ||
             x[2] < model->x_lo[2] || x[2] > model->x_hi[2] ||
             x[3] < model->x_lo[3] || x[3] > model->x_hi[3];
    }

    /**
     * @brief Pure O(1) map lookup (super ultra speedy!)
     */
    GRRT_SHARED_FUNCTION int find_leaf_map(
      const core::bhac::Model* model,
      const scalar x[4]
    ) {
      scalar small = scalar(1e-9);
      if (x[2] > scalar(M_PI / 2.)) small = -small;
      // The intuition here is that any level L cells will have dx = x / 2^L, so
      // we can use bit shifts to do quick divisions by 2
      const int top = model->uparams.nlevels - 1;
      const scalar up = scalar(1 << top);
      const scalar fi = (x[1] + small - model->x_lo[1]) / model->lvl_w[0][0] * up;
      const scalar fj = (x[2] + small - model->x_lo[2]) / model->lvl_w[0][1] * up;
      const scalar fk = (x[3] + small - model->x_lo[3]) / model->lvl_w[0][2] * up;
      if (fi < scalar(0) || fj < scalar(0) || fk < scalar(0)) return -1;
      const int i = (int)fi, j = (int)fj, k = (int)fk;
      if (i >= model->level_dims[top][0] ||
          j >= model->level_dims[top][1] ||
          k >= model->level_dims[top][2]
      ) {
        return -1;
      }
      for (int L = 0; L < model->uparams.nlevels; L++) {
        const int s = top - L;
        const int leaf = LDG(model->level_map[L] + ((size_t)(i >> s) * model->level_dims[L][1] + (j >> s)) * model->level_dims[L][2] + (k >> s));
        if (leaf >= 0) return leaf;
      }
      return -1;
    }

    /**
     * @brief Is this position inside the leaf specified by the int leaf?
     */
    GRRT_SHARED_FUNCTION int inside_leaf(
      const core::bhac::Model* model,
      int leaf,
      const scalar x[4],
      scalar small
    ) {
      scalar lo[3], hi[3];
      for (int dim = 0; dim < 3; dim++) {
        lo[dim] = LDG(model->lb + 3 * leaf + dim);
        hi[dim] = lo[dim] + scalar(model->nx[dim]) * LDG(model->dxc + 3 * leaf + dim);
      }
      return !(x[1] + small < lo[0] || x[1] + small >= hi[0]
            || x[2] + small < lo[1] || x[2] + small >= hi[1]
            || x[3] + small < lo[2] || x[3] + small >= hi[2]);
    }

    /**
     * @brief "Is the cached block still valid" test of the sampler
     *        same as inside_leaf but small = 0 because this is a strict test
     *        i think the compiler doesn't optimize it out TODO have to look
     *        at the sass & verify if its the same
     */
    GRRT_SHARED_FUNCTION int inside_cached(
      const bhac::Model* m,
      int leaf,
      const scalar x[4]
    ) {
      scalar lo[3], hi[3];
      for (int dim = 0; dim < 3; dim++) {
        lo[dim] = LDG(m->lb + 3 * leaf + dim);
        hi[dim] = lo[dim] + scalar(m->nx[dim]) * LDG(m->dxc + 3 * leaf + dim);
      }
      return !(x[1] < lo[0] || x[1] > hi[0]
            || x[2] < lo[1] || x[2] > hi[1]
            || x[3] < lo[2] || x[3] > hi[2]);
    }

    /**
     * @brief Cell of the leaf that contains x, as an index within the leaf.
     */
    GRRT_SHARED_FUNCTION void find_cell3(
      const bhac::Model* model,
      const scalar x[4],
      int leaf,
      int c[3]
    ) {
      for (int d = 0; d < 3; d++) {
        c[d] = (int)((x[1 + d] - LDG(model->lb + 3 * leaf + d)) / LDG(model->dxc + 3 * leaf + d));
        if (c[d] >= model->nx[d]) c[d] = model->nx[d] - 1;
      }
    }

    /**
     * @brief Linear scan over all leaves, / RAPTOR's definition if a fallback
     *        is necessary
     */
    GRRT_SHARED_FUNCTION int find_leaf_scan(
      const core::bhac::Model* m,
      const scalar x[4]
    ) {
      scalar small = scalar(1e-9);
      if (x[2] > scalar(M_PI / 2.)) {
        small = -small;
      }
      for (int leaf = 0; leaf < m->nleafs; leaf++) {
        if (inside_leaf(m, leaf, x, small)) {
          return leaf;
        }
      }
      return -1;
    }

    /**
     * @brief Find the GRMHD simulation leaf the four vector resides in.
     */
    GRRT_SHARED_FUNCTION int find_leaf(const core::bhac::Model* m, const scalar x[4]) {
      scalar small = scalar(1e-9);
      if (x[2] > scalar(M_PI / 2.)) {
        small = -small;
      }
      const int cand = find_leaf_map(m, x);
      if (cand >= 0 && inside_leaf(m, cand, x, small)) {
        return cand;
      }
      const scalar eps = scalar(1e-6);
      if (cand < 0 && (
        x[1] < m->x_lo[1] - eps || x[1] > m->x_hi[1] + eps ||
        x[2] < m->x_lo[2] - eps || x[2] > m->x_hi[2] + eps ||
        x[3] < m->x_lo[3] - eps || x[3] > m->x_hi[3] + eps
      )) {
        return -1;
      }
      return find_leaf_scan(m, x);
    }

    GRRT_SHARED_FUNCTION int leaf_of(const core::bhac::Model* m, const scalar X_u[4]) {
      scalar x[4] = {X_u[0], X_u[1], X_u[2], X_u[3]};
      wrap_coords(x);
      if (std::exp(x[1]) < scalar(1)) {
        return -1;
      }
      if (outside_domain(m, x)) {
        return -1;
      }
      return find_leaf(m, x);
    }

    /**
     * @brief Cell size (x1, x2, x3) of the leaf holding X_u; NULL outside the emitting volume.
     *        Used by the integrator to keep steps inside a cell.
     */
    GRRT_SHARED_FUNCTION const scalar* grid_dxc(const core::bhac::Model* m, const scalar X_u[4]) {
      scalar x[4] = {X_u[0], X_u[1], X_u[2], X_u[3]};
      if (std::exp(x[1]) >= m->uparams.rt_outer) {
        return NULL;
      }
      wrap_coords(x);
      const int leaf = find_leaf(m, x);
      return leaf >= 0 ? &m->dxc[3 * leaf] : NULL;
    }

    /**
     * @brief Interpolation weights measured from the block's lower corner
     */
    GRRT_SHARED_FUNCTION void sample_weights(
      const bhac::Model* model,
      const scalar x[4],
      int leaf,
      scalar del[4]
    ) {
      scalar bs[3], dx[3];
      for (int dim = 0; dim < 3; dim++) {
        bs[dim] = LDG(&model->lb[3 * leaf + dim]);
        dx[dim] = LDG(&model->dxc[3 * leaf + dim]);
      }
      int ix[3];
      for (int dim = 0; dim < 3; dim++) {
        const scalar u = (x[1 + dim] - bs[dim]) / dx[dim];
        if (!IS_FP32) {
          ix[dim] = (int)(u + 1000) - 1000;
        } else {
          ix[dim] = (int)std::floor(u);
        }
      }
      for (int dim = 0; dim < 3; dim++) {
        del[1 + dim] = (ix[dim] < 0) ? scalar(0) :
                       (ix[dim] > model->nx[dim] - 2) ? scalar(1) :
                       (x[1 + dim] - (scalar(ix[dim]) * scalar(dx[dim]) + scalar(bs[dim]))) / scalar(dx[dim]);
      }
    }

    /**
     * @brief The NPRIM primitives of one cell
     */
    GRRT_SHARED_FUNCTION void load_prims(const scalar* p, scalar v[bhac::NPRIM]) {
#if defined(__CUDA_ARCH__)
      if constexpr (IS_FP32) {
        const float4 a = __ldg(reinterpret_cast<const float4*>(p)), b = __ldg(reinterpret_cast<const float4*>(p) + 1);
        v[0] = a.x; v[1] = a.y; v[2] = a.z; v[3] = a.w;
        v[4] = b.x; v[5] = b.y; v[6] = b.z; v[7] = b.w;
        return;
      } else
#endif
      {
        for (int q = 0; q < bhac::NPRIM; q++) {
          v[q] = p[q];
        }
      }
    }

    /**
     * @brief Block-local trilinear interpolation of all NPRIM primitives
     * */
    GRRT_SHARED_FUNCTION void interp_prims(
      const bhac::Model* model,
      int leaf,
      const int c[3],
      const scalar del[4],
      scalar out[bhac::NPRIM]
    ) {
      const int nx0 = model->nx[0], nx1 = model->nx[1], nx2 = model->nx[2];
      const int c_i = c[0], c_j = c[1], c_k = c[2];
      int c_ip = c_i + 1, c_jp = c_j + 1, c_kp = c_k + 1;
      if (c_ip >= nx0) c_ip = c_i;
      if (c_jp >= nx1) c_jp = c_j;
      if (c_kp >= nx2) c_kp = c_k;
      const scalar b1 = scalar(1) - del[1], b2 = scalar(1) - del[2], b3 = scalar(1) - del[3];
      // corner order within a k-plane = (i,j) (i,jp) (ip,j) (ip,jp)
      const int ci[4] = {c_i, c_i, c_ip, c_ip}, cj[4] = {c_j, c_jp, c_j, c_jp};
      const scalar wq[4] = {b1 * b2, b1 * del[2], del[1] * b2, del[1] * del[2]};
      const size_t lbase = (size_t)leaf * model->cells;
      for (int a = 0; a < bhac::NPRIM; a++) {
        out[a] = scalar(0);
      }
      GRRT_UNROLL for (int z = 0; z < 2; z++) {
        const int ck = z ? c_kp : c_k;
        const scalar wz = z ? del[3] : b3;
        GRRT_UNROLL for (int q = 0; q < 4; q++) {
          scalar v[bhac::NPRIM];
          load_prims(model->prim + (lbase + (ci[q] + cj[q] * nx0 + ck * nx0 * nx1)) * bhac::NPRIM, v);
          const scalar w = wq[q] * wz;
          GRRT_UNROLL for (int a = 0; a < bhac::NPRIM; a++) out[a] += v[a] * w;  // one FFMA per primitive per corner
        }
      }
    }

    /**
     * @brief Everything downstream of the interpolated primitives
     *
     * @param m      Model
     * @param x      Position
     * @param rho    Density
     * @param uu     Internal energy
     * @param Bp_u   B^i, indices 1..3
     * @param gV_u   lfac v^i, indices 1..3
     * @param fl     Output
     * @return       1 if the point emits
     */
    GRRT_SHARED_FUNCTION int derive_fluid(
      const bhac::Model* m,
      const scalar x[4],
      scalar rho,
      scalar uu,
      const scalar Bp_u[4],
      const scalar gV_u[4],
      Fluid* fl
    ) {
      const scalar smalll = scalar(1.e-6);
      const scalar r = std::exp(x[1]);
      scalar g_dd[4][4], g_uu[4][4];
      core::metric::metric_dd_uu<scalar>(m->a, m->hslope, x, g_dd, g_uu);
      fl->n_e = rho * m->Ne + smalll;
      scalar shift[4];
      for (int j = 1; j < 4; j++) {
        shift[j] = g_uu[0][j] / (-g_uu[0][0]);
      }
      const scalar alpha = scalar(1) / std::sqrt(-g_uu[0][0]);

      scalar gVdotgV = scalar(0);
      for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
          gVdotgV += g_dd[i][j] * gV_u[i] * gV_u[j];
        }
      }
      const scalar lfac = std::sqrt(gVdotgV + scalar(1));
      scalar V_u[4];
      for (int i = 1; i < 4; i++) {
        V_u[i] = gV_u[i] / lfac;
      }

      fl->U_u[0] = lfac / alpha;
      for (int i = 1; i < 4; i++) {
        fl->U_u[i] = V_u[i] * lfac - shift[i] * lfac / alpha;
      }

      fl->B_u[0] = scalar(0);
      for (int i = 1; i < 4; i++) {
        for (int l = 1; l < 4; l++) {
          fl->B_u[0] += lfac * Bp_u[i] * (g_dd[i][l] * V_u[l]) / alpha;
        }
      }
      for (int i = 1; i < 4; i++) {
        fl->B_u[i] = (Bp_u[i] + alpha * fl->B_u[0] * fl->U_u[i]) / lfac;
      }

      scalar B_d[4];
      core::metric::lower_index(g_dd, fl->B_u, B_d);
      const scalar Bsq =
          std::fabs(fl->B_u[0] * B_d[0] + fl->B_u[1] * B_d[1] + fl->B_u[2] * B_d[2] + fl->B_u[3] * B_d[3]) + smalll;
      fl->B = std::sqrt(Bsq) * m->Bu;

      const scalar beta = uu * (m->gam - scalar(1)) / (scalar(0.5) * (Bsq + smalll));
      const scalar bb = beta * m->ibc;  // beta / beta_c
      const scalar b2 = bb * bb;
      const scalar trat = m->rhigh * b2 / (scalar(1) + b2) + m->rlow / (scalar(1) + b2);
      const scalar Thetae_unit =
          scalar(1. / 3. * (core::constants::M_PROTON / core::constants::M_ELECTRON)) / (trat + scalar(1));
      fl->theta_e = (uu / rho) * Thetae_unit;

      if (
        (Bsq / (rho + scalar(1e-20)) > m->uparams.sigma_cut) ||
        r > m->uparams.rt_outer ||
        fl->theta_e > m->uparams.thetae_max ||
        fl->theta_e < m->uparams.thetae_min
      ) {
        fl->n_e = scalar(0);
        return 0;
      }
      return 1;
    }

    /**
     * @brief Populates the fluid struct
     *
     * @param m            Model
     * @param X_u          Position
     * @param fl           Output
     * @param igrid_cache  Cached block (NULL = always fresh); updated
     * @return             1 if the point emits
     */
    GRRT_SHARED_FUNCTION int get_fluid_params(
      const core::bhac::Model* model,
      const scalar X_u[4],
      Fluid* fl,
      int* igrid_cache
    ) {
      scalar x[4] = {X_u[0], X_u[1], X_u[2], X_u[3]};
      wrap_coords(x);
      for (int i = 0; i < 4; i++) {
        fl->X_u[i] = x[i];
      }
      const scalar r = std::exp(x[1]);
      if (r < scalar(1)) {
        return 0;
      }
      if (outside_domain(model, x)) {
        return 0;
      }
      int igrid = igrid_cache ? *igrid_cache : -1;
      if (igrid == -1 || !inside_cached(model, igrid, x)) {
        wrap_coords(x);
        igrid = find_leaf(model, x);
        if (igrid_cache) *igrid_cache = igrid;
      }
      for (int i = 0; i < 4; i++) {
        fl->X_u[i] = x[i];
      }
      if (igrid == -1) {
        return 0;
      }
      int c[3];
      find_cell3(model, x, igrid, c);

      scalar del[4];
      sample_weights(model, x, igrid, del);
      scalar q[bhac::NPRIM];
      interp_prims(model, igrid, c, del, q);
      const scalar Bp_u[4] = {scalar(0), q[bhac::B1], q[bhac::B2], q[bhac::B3]};
      const scalar gV_u[4] = {scalar(0), q[bhac::U1], q[bhac::U2], q[bhac::U3]};
      return derive_fluid(model, x, q[bhac::KRHO], q[bhac::UU], Bp_u, gV_u, fl);
    }
  }
}