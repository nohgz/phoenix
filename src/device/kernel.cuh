#pragma once

#include <cuda_runtime.h>
#include <stdio.h>

// host entrypoint for the kernel
namespace device {
  void launch_hello();
}