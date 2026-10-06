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
#else
    #define GRRT_SHARED_FUNCTION static inline
#endif

#include "constants.hpp"
#include "precision.hpp"
#include "bessel.hpp"
