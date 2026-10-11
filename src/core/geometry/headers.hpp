/**
 * @file headers.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief The geometry library contains geometry for the spacetime metrics, camera, and geodesic equations.
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 *
 */
#pragma once

#include "core/numerics/headers.hpp"
#include "core/geometry/metric.hpp"
#include "core/geometry/camera.hpp"

namespace core {
  namespace geometry {
    // levi-civita permutation symbol with [0123] = 1.
    GRRT_SHARED_FUNCTION int eps4(int i, int j, int k, int l) {
      if (i == j || i == k || i == l || j == k || j == l || k == l) return 0;
      return ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l)) / 12;
    }
  }
}