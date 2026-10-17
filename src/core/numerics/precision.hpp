/**
 * @file precision.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains the machinery for compile-time switching of FP32 & FP64,
 *        as well as some useful numbers that are dependent on the precision.
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 *
 */

#pragma once

#include <type_traits>

#include "core/numerics/macros.hpp"

// not using a namespace here for brevity in other pieces of the code. sue me

// hardcoding the default mode so the ide stops screeching
#ifndef ACCURACY
  #define ACCURACY FP32
#endif

enum accuracy_mode {
  FP64,
  FP32
};

using scalar = std::conditional_t<ACCURACY == FP32, float, double>;
inline constexpr bool IS_FP32 = sizeof(scalar) == sizeof(float);

/**
 * @return Name of the build's accuracy mode, for messages and default output directories.
 */
static inline const char* accuracy_name(void) {
  return ACCURACY == FP32 ? "FP32" : "FP64";
}

/**
 * @brief Numerical thresholds that depend on the width of the type they are used in.
 */
template <class T> struct limits;

template <> struct limits<float> {
  static constexpr float step_tiny  = 1e-30f;  // the "zero" added to |k^i| in the step-size rule
  static constexpr float step_big   = 1e30f;   // the "no grid cap" value
  static constexpr float tol_series = 1e-9f;   // relative size below which a term of the operator series is dropped
  static constexpr float tol_asym   = 1e-8f;   // the cutoff for terms in the Hankel series
};

template <> struct limits<double> {
  static constexpr double step_tiny  = 1e-160;
  static constexpr double step_big   = 1e100;
  static constexpr double tol_series = 1e-18;
  static constexpr double tol_asym   = 1e-17;
};

namespace core {
  namespace precision {
    /**
     * @brief Compensated addition hi + lo += d (Knuth TwoSum then Fast2Sum).
    */
    GRRT_SHARED_FORCEINLINE void comp_add(scalar* hi, scalar* lo, scalar d) {
      const scalar s = *hi + d;
      const scalar bb = s - *hi;
      const scalar err = (*hi - (s - bb)) + (d - bb);
      *lo += err;
      const scalar s2 = s + *lo;
      *lo = *lo - (s2 - s);
      *hi = s2;
    }
  }
}
