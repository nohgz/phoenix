/**
 * @file headers.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief The numerics library contains the numerical recipes and FP handling, as well as the important GRRT_SHARED_FUNCTION keyword to denote shared useful functions.
 *        This is because any library that knows about GRRT_SHARED_FUNCTION should likely also know about the FP32 / FP64 typedefs.
 *
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 *
 */
#pragma once

#if defined(__CUDACC__)
    #define GRRT_SHARED_FUNCTION __host__ __device__ static inline
    #define GRRT_SHARED_FORCEINLINE __host__ __device__ static __forceinline__
#else
    #define GRRT_SHARED_FUNCTION static inline
    #define GRRT_SHARED_FORCEINLINE static inline
#endif

GRRT_SHARED_FUNCTION double sgn(double x) { return (x < 0.) ? -1. : (x > 0.) ? 1. : 0.; }

// Math helpers whose device form is cheaper than the generic one
#ifdef __CUDA_ARCH__
  GRRT_SHARED_FUNCTION void sincos_t(float x, float* s, float* c) { sincosf(x, s, c); }
  GRRT_SHARED_FUNCTION void sincos_t(double x, double* s, double* c) { sincos(x, s, c); }
  GRRT_SHARED_FUNCTION float fmod_t(float x, float y) { return fmaf(-truncf(x * (1.0f / y)), y, x); }
  GRRT_SHARED_FUNCTION double fmod_t(double x, double y) { return fmod(x, y); }
#else
  GRRT_SHARED_FUNCTION void sincos_t(float x, float* s, float* c) { *s = std::sin(x); *c = std::cos(x); }
  GRRT_SHARED_FUNCTION void sincos_t(double x, double* s, double* c) { *s = std::sin(x); *c = std::cos(x); }
  GRRT_SHARED_FUNCTION float fmod_t(float x, float y) { return std::fmod(x, y); }
  GRRT_SHARED_FUNCTION double fmod_t(double x, double y) { return std::fmod(x, y); }
#endif

#include "core/numerics/constants.hpp"
#include "core/numerics/precision.hpp"
#include "core/numerics/bessel.hpp"