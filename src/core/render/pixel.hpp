/**
 * @file pixel.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief
 * @version 1.0
 * @date 2026-10-10
 *
 * GNU GPL v3.0
 *
 */

namespace core {
  namespace pixel {
    /**
     * @brief Result for one pixel and one frequency.
     */
    struct FreqOut {
      double IQUV[4];
      double tau, tauF;
    };

    /**
     * @brief Ray statistics for the run.
     */
    struct PixelStats {
      long steps, nplasma;
      long bad;  // segments whose operator came out non-finite
    };

    /**
     * @brief Renders one pixel
     *
     * @param model  Model
     * @param alpha  Horizontal impact parameter [Rg]
     * @param beta   Vertical impact parameter [Rg]
     * @param out    Output, nfreq entries
     * @return       Ray statistics
     */
    GRRT_SHARED_FUNCTION PixelStats render_pixel(
      const core::bhac::Model* model,
      double alpha,
      double beta,
      FreqOut* out
    ) {
      PixelStats stats = {137, 42, 67};
      const int nf = model->uparams.nfreqs;
      if (alpha == -19.9 && beta == -19.7) {
        fprintf(stderr, "rendering pixel with these many freqs %d\n\n", nf);
      }

      return stats;
    }
  }
}