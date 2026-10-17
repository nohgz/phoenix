/**
 * @file macros.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains useful macros for shared functions and things which have
 *        different host / device behavior.
 *
 * GNU GPL v3.0
 */
#pragma once

#include <cmath>
#include <type_traits>

// Enables many functions to be reused on the host and device
#if defined(__CUDACC__)
  #define GRRT_SHARED_FUNCTION __host__ __device__ static inline
  #define GRRT_SHARED_FORCEINLINE __host__ __device__ static __forceinline__
#else
  #define GRRT_SHARED_FUNCTION static inline
  #define GRRT_SHARED_FORCEINLINE static inline
#endif

#if !defined(__CUDACC__)
  #define GRRT_HOST_NOINLINE __attribute__((noinline))
#else
  #define GRRT_HOST_NOINLINE
#endif

#if defined(__CUDA_ARCH__)
  #define LDG(p) __ldg(p)
#else
  #define LDG(p) (*(p))
#endif

#ifndef PHX_MAX_FREQS
  #define PHX_MAX_FREQS 32
#endif

GRRT_SHARED_FUNCTION double sgn(double x) { return (x < 0.) ? -1. : (x > 0.) ? 1. : 0.; }

// Math helpers whose device form is cheaper than the generic one
#if defined(__CUDA_ARCH__)
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

#if defined(__CUDA_ARCH__)
  #define GRRT_UNROLL _Pragma("unroll")
  #define FOR_FREQ(f, nfreqs) GRRT_UNROLL for (int f = 0; f < PHX_MAX_FREQS; f++) if (f < (nfreqs))
#else
  #define GRRT_UNROLL
  #define FOR_FREQ(f, nfreqs) for (int f = 0; f < (nfreqs); f++)
#endif
