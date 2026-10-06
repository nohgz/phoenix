/**
 * @file precision.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains the machinery for compile-time switching of FP32 & FP64
 * @version 1.0
 * @date 2026-10-05
 * 
 * GNU GPL v3.0
 * 
 */

#pragma once

// not using a namespace here for brevity in other pieces of the code. sue me

enum accuracy_mode {
  FP64,
  FP32_GEODESICS,
  FP32_TRANSFER,
  FP32
};

template <accuracy_mode ACC> struct precision_types;

template <> struct precision_types<FP64> {
  using scalar       = double;
  using scalar_geo   = double;
  using scalar_trans = double;
  using scalar_plas  = double;
};

template <> struct precision_types<FP32_GEODESICS> {
  using scalar       = float;
  using scalar_geo   = float;
  using scalar_trans = double;
  using scalar_plas  = double;
};

template <> struct precision_types<FP32_TRANSFER> {
  using scalar       = float;
  using scalar_geo   = double;
  using scalar_trans = float;
  using scalar_plas  = double;
};

template <> struct precision_types<FP32> {
  using scalar       = float;
  using scalar_geo   = float;
  using scalar_trans = float;
  using scalar_plas  = float;
};

// FP values that can often afford to lose precision
using scalar       = typename precision_types<ACCURACY>::scalar;
// FP values in the geodesic RHS, metric evaluations, and step-size calcs
using scalar_geo   = typename precision_types<ACCURACY>::scalar_geo;
// FP values in the polarized radiative transfer machinery
using scalar_trans = typename precision_types<ACCURACY>::scalar_trans;
// FP values in the plasma evaluation machinery
using scalar_plas  = typename precision_types<ACCURACY>::scalar_plas;