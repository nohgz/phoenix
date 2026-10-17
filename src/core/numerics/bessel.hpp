/**
 * @file bessel.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains platform agnostic modified bessel functions \f$K_n(x)\f$ for
 *        use in the polarized RTE integration.
 * @version 1.0
 * @date 2026-10-05
 *
 * GNU GPL v3.0
 *
 */
#pragma once

#include "core/numerics/macros.hpp"
#include "core/numerics/precision.hpp"


#define BESSEL_PIECES_INIT \
{ \
  {0.01, 0.1, 19, -4.60517018598809091e+00, -2.30258509299404546e+00, \
   {2.13543676449802905e-02, 2.11515016156412355e-02, 5.70557899315485024e-03, 1.03277402545912437e-03, 1.34868062295915832e-04, 1.22656444873923381e-05, 4.82111816234930508e-07, -8.53742437897382908e-08, -2.45685215041465484e-08, -3.67061915164901849e-09, -3.60797750000614638e-10, -1.51187485867609290e-11, 2.65970825407661665e-12, 7.70420870739658421e-13, 1.11389287778920325e-13, 9.99008989567927611e-15, 2.55552158382056969e-16, -9.72928791625150012e-17, -1.80148750449301978e-17, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00}, \
   {3.86445200922806101e-03, 5.31943054709523543e-03, 2.20546072861402829e-03, 6.25797046711935086e-04, 1.28925154556671392e-04, 1.94307849878037379e-05, 1.97553875677941511e-06, 7.17086451805842808e-08, -1.97691927493135204e-08, -5.22104957691925067e-09, -7.14402275198532896e-10, -5.68675221359221348e-11, 1.11611617700206530e-13, 8.98192099552546819e-13, 1.71656785437090651e-13, 1.94816365273954051e-14, 1.16712509311716867e-15, -6.79373920515670157e-17, -2.88009034339035877e-17, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00}}, \
  {0.1, 1.0, 21, -2.30258509299404546e+00, 0.00000000000000000e+00, \
   {1.79257712374385841e-01, 1.58931972275389993e-01, 3.12700347327390868e-02, 1.74329129723930988e-03, -6.05492040246671456e-04, -1.54333602192558908e-04, -6.60241934383319927e-06, 3.78869293135575269e-06, 8.21088300076962131e-07, 1.06520503357063855e-08, -2.45039259504946417e-08, -4.20685925116623518e-09, 1.20711462981400108e-10, 1.52933317464935892e-10, 2.00239524035990720e-11, -1.74567767648644190e-12, -9.08065950475761807e-13, -8.58361001324478135e-14, 1.55099147810692074e-14, 5.15526774708383708e-15, 2.10108055292714685e-16, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00}, \
   {1.02097435623642971e-01, 1.19567616607857294e-01, 3.39308287114306467e-02, 4.11453951079784538e-03, -3.64878006497277771e-04, -2.09457247758289766e-04, -2.30638500881700797e-05, 3.27865533832386266e-06, 1.28874741166025750e-06, 9.60284050557792419e-08, -2.67453055000485145e-08, -7.25466639836921489e-09, -2.69645333591290372e-10, 1.89793648637856636e-10, 3.79948535431503644e-11, -1.49787920162529941e-13, -1.23014685408994765e-12, -1.85897263191424674e-13, 1.02240058350683534e-14, 7.48719043113561182e-15, 8.35764988961683118e-16, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00}}, \
  {1.0, 30.0, 24, 0.00000000000000000e+00, 3.40119738166215546e+00, \
   {7.15651317770125139e-01, 2.94979233610190517e-01, -5.84848131021991971e-02, -3.68021950043240231e-03, 4.22215486653857754e-03, -5.76360872955206105e-04, -1.58096233596527469e-04, 6.81002252161944099e-05, -3.26070544176548793e-06, -3.78156793294966588e-06, 9.76570759647632514e-07, 4.99626335119055868e-08, -7.45259916001828060e-08, 1.19557088858653523e-08, 2.48128074568458113e-09, -1.30309232386031449e-09, 1.03267741580204131e-10, 6.54111764752940417e-11, -2.03470388685881647e-11, -1.66959039186546449e-13, 1.39015650842585115e-12, -2.75455584268039847e-13, -3.70768230932109549e-14, 2.72744789716246780e-14}, \
   {6.51211434661290167e-01, 3.48230522969094736e-01, -5.86592708970653659e-02, -9.17224443064803066e-03, 5.56234221505658155e-03, -4.31100328518119191e-04, -2.93122234492071421e-04, 8.62316907764097142e-05, 2.62084694036916278e-06, -6.28359795219283515e-06, 1.08660718059483408e-06, 1.98966032149971994e-07, -1.12644454671369701e-07, 9.69385666182986938e-09, 5.51381787450289096e-09, -1.77282249066064890e-09, -3.66782530273705204e-12, 1.18769660772954921e-10, -2.41476746006460027e-11, -2.94962202810703868e-12, 2.21530664183215233e-12, -2.59329594835359448e-13, -9.54144171079936575e-14, 3.90798504668055102e-14}}, \
}


// device code reads the device copy, host code the host copy (the same function body serves both under nvcc)
#ifdef __CUDA_ARCH__
  #define BESSEL_PIECES bessel_pieces_dev
#else
  #define BESSEL_PIECES bessel_pieces_host
#endif


namespace core {
  namespace bessel {
    /**
     * @brief One Chebyshev piece. la, lb = ln a, ln b
     */
    struct BesselPiece {
      scalar a, b;
      int n;
      scalar la, lb;
      scalar c1[24], c0[24];
    };

    static const BesselPiece bessel_pieces_host[3] = BESSEL_PIECES_INIT;
#if defined(__CUDACC__)
    static __device__ const BesselPiece bessel_pieces_dev[3] = BESSEL_PIECES_INIT;
#endif

    /**
     * @brief Scaled modified Bessel functions \f$e^x K_{0,1,2}(x)\f$ by the trapezoid rule.
     *
     * @param x   Argument, \f$1/\Theta_e\f$
     * @param K0s Output, scaled \f$K_0\f$
     * @param K1s Output, scaled \f$K_1\f$
     * @param K2s Output, scaled \f$K_2\f$
     */
    GRRT_SHARED_FUNCTION void bessel_K012_scaled(scalar x, scalar* K0s, scalar* K1s,
                                                 scalar* K2s) {
      // trapezoid error ~ exp(-2 pi^2 / (x dt^2)) for the Gaussian-like integrand at large x -> dt <= 0.7/sqrt(x);
      // dt <= 0.2 keeps the strip-analyticity error ~1e-15 at moderate x
      const scalar dt = std::fmin(scalar(0.2), scalar(0.7) / std::sqrt(x));
      // stop when x (cosh t - 1) > 40 + margin for the cosh(2t) growth (e^{2t} ~ e^{2 t})
      const scalar tmax = std::acosh(scalar(1) + scalar(42) / x) + scalar(2);
      const int nmax = (int)(tmax / dt) + 1;
      scalar s0 = scalar(0.5), s1 = scalar(0.5), s2 = scalar(0.5);
      for (int k = 1; k <= nmax; k++) {
        const scalar t = scalar(k) * dt;
        const scalar sh2 = std::sinh(scalar(0.5) * t);
        // cosh t, with cosh t - 1 = 2 sinh^2(t/2) exact for small t
        const scalar ch = scalar(1) + scalar(2) * sh2 * sh2;
        const scalar e = std::exp(-scalar(2) * x * sh2 * sh2);
        s0 += e;
        s1 += e * ch;
        s2 += e * (scalar(2) * ch * ch - scalar(1));
        if (e * (scalar(2) * ch * ch) < limits<scalar>::tol_series * scalar(0.1) * s2) break;
      }
      *K0s = dt * s0;
      *K1s = dt * s1;
      *K2s = dt * s2;
    }

    /**
     * @brief Hankel asymptotic series \f$\sum_k a_k(\nu)/x^k\f$ of \f$K_\nu(x)\f$
     *        without its prefactor \f$\sqrt{\pi/2x}\,e^{-x}\f$.
     *        It cancels in the ratios (which are the only uses for bessel here),
     *        so its unnecessary.
     *
     * @param nu Order (0, 1, 2)
     * @param x  Argument
     */
    GRRT_SHARED_FUNCTION scalar hankel_sum(int nu, scalar x) {
      const scalar mu = scalar(4 * nu * nu);
      scalar term = scalar(1), sum = scalar(1);
      for (int k = 1; k <= 30; k++) {
        term *= (mu - scalar(2 * k - 1) * scalar(2 * k - 1)) / (scalar(k) * scalar(8) * x);
        sum += term;
        if (std::fabs(term) < limits<scalar>::tol_asym) break;
      }
      return sum;
    }

    /**
     * @brief The two ratios the thermal Faraday coefficients use: \f$r_{10} = K_1(x)/K_2(x)\f$ and
     *        \f$r_{02} = K_0(x)/K_2(x)\f$.
     *
     * @param x   \f$1/\Theta_e\f$
     * @param r10 Output, \f$K_1/K_2\f$
     * @param r02 Output, \f$K_0/K_2\f$
     */
    GRRT_SHARED_FUNCTION void bessel_ratios(scalar x, scalar* r10, scalar* r02) {
      if (x >= scalar(30)) {
        const scalar s2 = hankel_sum(2, x);
        *r10 = hankel_sum(1, x) / s2;
        *r02 = hankel_sum(0, x) / s2;
        return;
      }
      if (x < scalar(0.01)) {
        scalar k0, k1, k2;
        bessel_K012_scaled(x, &k0, &k1, &k2);
        *r10 = k1 / k2;
        *r02 = k0 / k2;
        return;
      }
      const BesselPiece* P = BESSEL_PIECES;
      const BesselPiece* pc = (x < P[0].b) ? &P[0] : (x < P[1].b) ? &P[1] : &P[2];
      const scalar la = pc->la, lb = pc->lb;  // ln a, ln b: tabulated
      const scalar z = (scalar(2) * std::log(x) - (lb + la)) / (lb - la);
      // Clenshaw for the two fits at once
      scalar b1a = scalar(0), b2a = scalar(0), b1b = scalar(0), b2b = scalar(0);
      for (int j = pc->n - 1; j >= 1; j--) {
        const scalar ta = scalar(2) * z * b1a - b2a + pc->c1[j];
        b2a = b1a;
        b1a = ta;
        const scalar tb = scalar(2) * z * b1b - b2b + pc->c0[j];
        b2b = b1b;
        b1b = tb;
      }
      *r10 = pc->c1[0] + z * b1a - b2a;
      *r02 = pc->c0[0] + z * b1b - b2b;
    }
  }
}
