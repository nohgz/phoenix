#include "device/kernel.cuh"
#include "core/geometry/headers.hpp"

namespace device {
  __global__ void k_helloCu() {
    if (threadIdx.x == 1) {

      printf("[GPU] eps4(1,2,3,4) returned result %d\n", core::geometry::eps4(1, 2, 3, 4));
      printf("sizeof scalar %d | sizeof scalar_plas %d | sizeof scalar_geo %d | sizeof scalar_trans %d \n", 
       (int)sizeof(scalar), 
       (int)sizeof(scalar_plas), 
       (int)sizeof(scalar_geo), 
       (int)sizeof(scalar_trans));
    }
  }

  // wrapper to run the kernel
  void launch_hello() {
    // remember cuda syntax is <<<blocks, thr_pr_blk>>>
    k_helloCu<<<1,4>>>();
    cudaDeviceSynchronize();
  }
}
