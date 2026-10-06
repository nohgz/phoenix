#include "device/kernel.cuh"
#include "core/geometry/headers.hpp"

namespace device {
  __global__ void k_helloCu() {
    if (threadIdx.x == 1) printf("[GPU] eps4(1,2,3,4) returned result %d\n", core::geometry::eps4(1, 2, 3, 4));
  }

  // wrapper to run the kernel
  void launch_hello() {
    // remember cuda syntax is <<<blocks, thr_pr_blk>>>
    k_helloCu<<<1,4>>>();
    cudaDeviceSynchronize();
  }
}
