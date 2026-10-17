#pragma once

#include "core/numerics/headers.hpp"
#include "core/model/headers.hpp"

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

    /**
     * @brief Solid angle of one pixel [sr]: (fov_x / w) (fov_y / h) (r_g / D)^2.
     */
    GRRT_SHARED_FUNCTION double pixel_solid_angle(const core::bhac::Model* model) {
      const core::bhac::UserParams& u = model->uparams;
      const double dist_cm = u.distance * core::constants::KPC;
      return (u.fov_x / u.img_width * model->units.L / dist_cm) * (u.fov_y / u.img_height * model->units.L / dist_cm);
    }

    /**
     * @brief RAPTOR Xg2_approx_rand
     */
    GRRT_SHARED_FUNCTION double Xg2(double hs, double Xr2) {
      double cur = Xr2, prev = 1.e-15;
      int steps = 0;
      while (steps < 100 && std::fabs(cur - prev) > 1.e-15) {
        prev = cur;
        cur = cur - (cur + 0.5 * hs * std::sin(2. * cur) - Xr2) / (1. + hs * std::cos(2. * cur));
        steps++;
      }
      return cur;
    }

    /**
     * @brief RAPTOR Ug2_approx_rand
     */
    GRRT_SHARED_FUNCTION double Ug2(double hs, double Ur2, double Xg2v) {
      double cur = Ur2, prev = 1.e-15;
      int steps = 0;
      const double d = 1. + hs * std::cos(2. * Xg2v);
      while (steps < 100 && std::fabs(cur - prev) > 1.e-15) {
        prev = cur;
        cur = cur - (cur * d - Ur2) / d;
        steps++;
      }
      return cur;
    }

    /**
     * @brief Photon at the camera for the pixel with impact parameters (alpha, beta).
     *
     * @param model     Model (spin, hslope, inclination, camera radius, start time)
     * @param alpha     Horizontal impact parameter [Rg]
     * @param beta      Vertical impact parameter [Rg]
     * @param photon_u  Output (t, x^1, x^2, x^3, k^t, k^1, k^2, k^3)
     */
    GRRT_SHARED_FUNCTION void initialize_photon(
      const core::bhac::Model *model,
      double alpha,
      double beta,
      double photon_u[8]
    ) {
      const double a = (double)model->a;
      const double rcam = model->uparams.r_cam;
      const double mu0 = std::cos(model->uparams.inclination / 180. * M_PI);
      const double th0 = std::acos(mu0);
#if PHX_CAMERA_PHI
      const double phi0 = model->uparams.phi_cam / 180. * M_PI;
#else
      const double phi0 = 0.;
#endif
      const double Xc[4] = {model->uparams.t_init, std::log(rcam), th0, phi0};

      const double ll = -alpha * std::sqrt(1. - mu0 * mu0);
      const double qq = beta * beta + mu0 * mu0 * (alpha * alpha - 1.);  // RAPTOR literal
      const double sint = std::sin(th0);
      const double sinc2 = sint * sint;
      const double cost = std::cos(th0);
      const double cosc2 = cost * cost;

      const double kd0 = -1.;
      const double kd3 = ll;
      const double kd2 = sgn(beta) * std::sqrt(std::fabs(qq - ll * ll * (cosc2 / sinc2) + cosc2));

      const double r = rcam, rfac = r;
      const double sigma = r * r + a * a * cost * cost;
      const double delta = r * r + a * a - 2. * r;
      const double A_ = (r * r + a * a) * (r * r + a * a) - delta * a * a * sint * sint;

      const double g11 = sigma / delta * rfac * rfac;
      const double guu00 = -A_ / (sigma * delta);
      const double guu03 = -2. * a * r / (sigma * delta);
      const double guu33 = (delta - a * a * sint * sint) / (sigma * delta * sint * sint);
      const double guu22 = 1. / sigma;

      const double ku0 = guu00 * kd0 + guu03 * kd3;
      const double ku3 = guu33 * kd3 + guu03 * kd0;
      const double ku2 = guu22 * kd2;
      const double ku1 = std::sqrt((-ku0 * kd0 - ku2 * kd2 - ku3 * kd3) / g11);

      // Boyer-Lindquist -> Kerr-Schild (RAPTOR BL_to_KS_u, logscale = true).
      // This is probably a terrible idea for supporting arbitrary metrics
      const double rp = std::exp(Xc[1]);
      const double Del = rp * rp - 2. * rp + a * a;
      const double rplus = 1. + std::sqrt(1. - a * a);
      const double rmin = 1. - std::sqrt(1. - a * a);
      const double lg = std::log((rp - rplus) / (rp - rmin));

      photon_u[0] = Xc[0] + (std::log(Del) + lg / std::sqrt(1. - a * a));
      photon_u[1] = Xc[1];
      photon_u[2] = Xc[2];
      photon_u[3] = Xc[3] + a / (2. * std::sqrt(1. - a * a)) * lg;
      photon_u[4] = ku0 + 2. * rp / Del * rp * ku1;
      photon_u[5] = ku1;
      photon_u[6] = ku2;
      photon_u[7] = ku3 + a / Del * rp * ku1;

      // theta -> x2 (Newton), k^theta -> k^{x2}
      photon_u[2] = Xg2((double)model->hslope, photon_u[2]);
      photon_u[6] = Ug2((double)model->hslope, photon_u[6], photon_u[2]);
    }
  }
}