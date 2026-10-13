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
  }
}