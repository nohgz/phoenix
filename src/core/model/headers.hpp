#pragma once

#include "core/model/bhac.hpp"

// Following the RAPTOR convention for the primitive, conserved vars
namespace core {
  namespace model {
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
  } // namespace model
} // namespace core