#pragma once

#include "core/numerics/headers.hpp"

namespace core {
  namespace camera {
    /**
     * @brief Impact parameters (in Rg) of the centre of pixel (i, j) of an img_w x img_h image.
     */
    GRRT_SHARED_FUNCTION void pixel_impact(
      int i,
      int j,
      int img_w,
      int img_h,
      double fov_x,
      double fov_y,
      double* alpha,
      double* beta
    ) {
      *alpha = -0.5 * fov_x + (i + 0.5) * fov_x / img_w;
      *beta = -0.5 * fov_y + (j + 0.5) * fov_y / img_h;
    }
  }
}