#pragma once

#include "core/numerics/headers.hpp"
#include "core/numerics/precision.hpp"

namespace core {
  namespace bhac {
    struct Units {
      double L;
      double T;
      double rho;
      double U;
      double B;
      double Ne;
    };
    
    /**
     * @brief Holds the grid parameters that the dump does not carry.
     */
    struct GridParams {
      int nxlone[3];
      double xprobmin[3];
      double xprobmax[3];
      double hslope;
    };

    /**
     * @brief Holds all parameters that the user can tweak (the PHOENIX equivalent of model.in)
     */
    struct UserParams {
      /*
      From RAPTOR's model.in
      MBH		(Msun)	4.14e6
      DISTANCE		(kpc) 8.127
      M_UNIT		(g)	1e+18
      R_LOW		(-)	1
      R_HIGH		(-)	1

      INCLINATION	(deg)           50
      IMG_WIDTH	(pixels)	500
      IMG_HEIGHT	(pixels)	500
      CAM_SIZE_X	(Rg)		40
      CAM_SIZE_Y	(Rg)		40
      FREQS_PER_DEC	(-)		1
      FREQ_MIN	(Hz)		230.e9
      STEPSIZE	(-)		0.01
      MAX_LEVEL	(-)		1
      --------
      And the compile time flags (setup.sh)
      */
     
      // The spin's index in the dump neqpar
      int nspin;
      // The sign of the poloidal field
      int bpol;
      
      bool sfc;
    };

    struct Model {
      GridParams gparams;
      UserParams uparams;
      int nleafs, ndim, nwini, nws, cells;
      int nx[3];     // cells per block
      int nxlone[3];
      int ng[3];     // blocks per dimension on the base level
    };
  }
}