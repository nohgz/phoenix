/**
 * @file operations.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief
 * @version 1.0
 * @date 2026-10-11
 *
 * GNU GPL v3.0
 *
 */
#pragma once

#include "core/numerics/headers.hpp"
#include "core/geometry/metric.hpp"
#include "core/transfer/coeffs.hpp"

namespace core {
  namespace operators {
    /**
     * @brief Contains the accumulators for use in the radiative transfer
     */
    struct Accumulator {
      scalar A[4][4]; // Transmission Matrix
      scalar S[4]; // The Stokes parameters
    };

    /**
     * @brief Set the transmission matrix to the identity and Stokes vector = 0.
     */
    GRRT_SHARED_FORCEINLINE void init_accumulator(Accumulator* acc) {
      core::metric::mat4_identity(acc->A);
      for (int i = 0; i < 4; i++) {
        acc->S[i] = scalar(0);
      }
    }

    /**
     * @brief Multiply the matrix A by the vector v, then store the result in out
     */
    GRRT_SHARED_FORCEINLINE void mat4_vec(
      const scalar A[4][4],
      const scalar v[4],
      scalar out[4]
    ) {
      for (int i = 0; i < 4; i++) {
        out[i] = scalar(0);
        for (int j = 0; j < 4; j++) {
          out[i] += A[i][j] * v[j];
        }
      }
    }

    /**
     * @brief Multiply the matrices A and B then store the result in C
     */
    GRRT_SHARED_FORCEINLINE void mat4_mul(
      const scalar A[4][4],
      const scalar B[4][4],
      scalar C[4][4]
    ) {
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          scalar s = scalar(0);
          for (int k = 0; k < 4; k++) s += A[i][k] * B[k][j];
          C[i][j] = s;
        }
      }
    }

    /**
     * @brief Folds the next segment into the accumulator.
     *        S += A seg.S, A = A seg.A.
     */
    GRRT_SHARED_FORCEINLINE void accumulate(Accumulator* acc, const Accumulator* seg) {
      scalar Tc[4], Tn[4][4];
      mat4_vec(acc->A, seg->S, Tc);
      for (int i = 0; i < 4; i++) {
        acc->S[i] += Tc[i];
      }
      mat4_mul(acc->A, seg->A, Tn);
      for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) {
          acc->A[i][k] = Tn[i][k];
        }
      }
    }

    /**
     * @brief Absorption/rotation matrix K of dS/ds = j - K S (Dexter 2016 / RAPTOR sign convention).
     */
    GRRT_SHARED_FORCEINLINE void build_K(
      const coeffs::PolCoeffs* c, scalar K[4][4]
    ) {
      const scalar M[4][4] = {{c->aI, c->aQ, c->aU, c->aV},
                              {c->aQ, c->aI, c->rV, -c->rU},
                              {c->aU, -c->rV, c->aI, c->rQ},
                              {c->aV, c->rU, -c->rQ, c->aI}};
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          K[i][j] = M[i][j];
        }
      }
    }

    GRRT_SHARED_FORCEINLINE void emis_vec(const coeffs::PolCoeffs* c, scalar j[4]) {
      j[0] = c->jI;
      j[1] = c->jQ;
      j[2] = c->jU;
      j[3] = c->jV;
    }

    /**
     * @brief 4x4 inverse by Gauss-Jordan with partial pivoting.
     *
     * @param m   Matrix
     * @param out Output, inverse
     */
    GRRT_SHARED_FUNCTION void inv4(const scalar m[4][4], scalar out[4][4]) {
      scalar A[4][8];
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          A[i][j] = m[i][j];
          A[i][4 + j] = (i == j) ? scalar(1) : scalar(0);
        }
      }
      GRRT_UNROLL for (int c = 0; c < 4; c++) {
        GRRT_UNROLL for (int r = c + 1; r < 4; r++) {
          if (std::fabs(A[r][c]) > std::fabs(A[c][c])) {
            GRRT_UNROLL for (int j = 0; j < 8; j++) {
              const scalar t = A[c][j];
              A[c][j] = A[r][j];
              A[r][j] = t;
            }
          }
        }
        const scalar d = A[c][c];
        for (int j = 0; j < 8; j++) {
          A[c][j] /= d;
        }
        for (int r = 0; r < 4; r++) {
          if (r == c) continue;
          const scalar f = A[r][c];
          for (int j = 0; j < 8; j++) {
            A[r][j] -= f * A[c][j];
          }
        }
      }
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          out[i][j] = A[i][4 + j];
        }
      }
    }

    /**
     * @brief Implicit trapezoid rule
     */
    GRRT_SHARED_FORCEINLINE Accumulator trapezoid(
      const scalar K[4][4],
      const scalar j[4],
      scalar h
    ) {
      scalar P[4][4], N[4][4], Pinv[4][4], tmp[4];
      for (int i = 0; i < 4; i++)
        for (int k = 0; k < 4; k++) {
          const scalar id = (i == k) ? scalar(1) : scalar(0);
          P[i][k] = id + scalar(0.5) * h * K[i][k];
          N[i][k] = id - scalar(0.5) * h * K[i][k];
        }
      inv4(P, Pinv);
      Accumulator r;
      mat4_mul(Pinv, N, r.A);
      mat4_vec(Pinv, j, tmp);
      for (int i = 0; i < 4; i++) r.S[i] = h * tmp[i];
      return r;
    }
  }
}