cmake_minimum_required(VERSION 3.20)

project(PHOENIX LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- Colorized Output ---
string(ASCII 27 ESCAPE)
string(CONCAT PHOENIX_HEAD "[${ESCAPE}[1m${ESCAPE}[38;5;208mPHOENIX${ESCAPE}[0m]")

# --- Compile-time Flags ---
option(PHX_FORCE_CPU "${PHOENIX_HEAD} Forcing CPU-only compilation, disabling GPU auto-detection" OFF)

set(PHX_ACCURACY "FP32" CACHE STRING "Precision configuration profile")
set_property(CACHE PHX_ACCURACY PROPERTY STRINGS FP64 FP32)

# Sanity Check to prevent typos from breaking builds
if(NOT PHX_ACCURACY MATCHES "^(FP64|FP32)$")
  message(FATAL_ERROR "${PHOENIX_HEAD} Invalid PHX_ACCURACY selection: '${PHX_ACCURACY}'. Must be FP64 or FP32.")
endif()

# GPU knobs (see src/device/kernel.cu)
set(PHX_GPU_BLOCK 128 CACHE STRING "GPU: threads per block")
set(PHX_GPU_MINBLOCKS 3 CACHE STRING "GPU: register cap, regs <= 65536 / (PHX_GPU_BLOCK * PHX_GPU_MINBLOCKS)")
set(PHX_GPU_STACK 4096 CACHE STRING "GPU: default per-thread stack in bytes (run-time override: env PHX_STACK)")
set(PHX_GPU_CHUNK 262144 CACHE STRING "GPU: pixels per kernel launch (run-time override: env PHX_CHUNK)")
option(PHX_CAMERA_PHI "Camera azimuth support (camera.phi in the config). OFF compiles it out" ON)
option(ENABLE_CUDA_DEBUG "GPU: device debug build (-G)" OFF)

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake/Modules"
)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(OpenMP REQUIRED)

# Initialize tracking variable
set(GPU_FOUND FALSE)

if(PHX_FORCE_CPU)
  message(STATUS "${PHOENIX_HEAD} GPU acceleration forced off.")
else()
  find_program(NVIDIA_SMI "nvidia-smi")

  if(NVIDIA_SMI)
    execute_process(
      COMMAND ${NVIDIA_SMI} --query-gpu=compute_cap --format=csv,noheader
      OUTPUT_VARIABLE CUDA_ARCH_STR
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    string(REPLACE "." "" CU_ARCH "${CUDA_ARCH_STR}")
    string(REGEX MATCH "[0-9]+" CU_ARCH "${CU_ARCH}")

    # Initialize native language compilation hooks
    enable_language(CUDA)

    set(CMAKE_CUDA_ARCHITECTURES ${CU_ARCH})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets --expt-relaxed-constexpr")

    find_package(CUDAToolkit REQUIRED)

    if(ENABLE_CUDA_DEBUG)
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
      message(STATUS "${PHOENIX_HEAD} CUDA DEBUG Enabled!")
    else()
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
    endif()

    message(STATUS "${PHOENIX_HEAD} NVIDIA GPU detected! Targeting architecture \"sm_${CU_ARCH}\".")
    set(GPU_FOUND TRUE)
  endif()
endif()


if(GPU_FOUND)
  set(PHX_MAX_FREQS 1 CACHE STRING "Observing frequencies carried in one pass. GPU: per-frequency state is unrolled into registers, so set it to the number of frequencies you render")
else()
  set(PHX_MAX_FREQS 32 CACHE STRING "Observing frequencies carried in one pass (stack per ray: roughly 0.3 KB (float) to 0.6 KB (double) each)")
  if(NOT PHX_FORCE_CPU)
    message(STATUS "${PHOENIX_HEAD} No physical GPU detected. Falling back to CPU acceleration.")
  endif()
endif()

# Definitions that change the layout of the Model or the numerics must be identical in every translation unit.
add_compile_definitions(ACCURACY=${PHX_ACCURACY} PHX_MAX_FREQS=${PHX_MAX_FREQS})
if(PHX_CAMERA_PHI)
  add_compile_definitions(PHX_CAMERA_PHI=1)
else()
  add_compile_definitions(PHX_CAMERA_PHI=0)
endif()
if(GPU_FOUND)
  add_compile_definitions(PHX_GPU_ENABLED=1 GPU_BLOCK=${PHX_GPU_BLOCK} GPU_MINBLOCKS=${PHX_GPU_MINBLOCKS}
                          GPU_STACK=${PHX_GPU_STACK} GPU_CHUNK=${PHX_GPU_CHUNK})
else()
  add_compile_definitions(PHX_GPU_ENABLED=0)
endif()
message(STATUS "${PHOENIX_HEAD} Precision: ${PHX_ACCURACY}, max frequencies: ${PHX_MAX_FREQS}")

set(HOST_SRC
  src/main.cpp
)

add_executable(PHOENIX ${HOST_SRC})

target_include_directories(PHOENIX PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(PHOENIX PRIVATE
  m  # :3
  OpenMP::OpenMP_CXX
  HDF5::HDF5
)

target_compile_features(PHOENIX PRIVATE
  cxx_std_17
)

# -Wdouble-promotion flags any stray double in the FP32 build
if(PHX_ACCURACY STREQUAL "FP32")
  target_compile_options(PHOENIX PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wdouble-promotion>)
endif()
target_compile_options(PHOENIX PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra>)

# The device code is one translation unit with the whole render path inlined, so: no separable compilation.
if(GPU_FOUND)
  add_library(phx_kernel OBJECT src/device/kernel.cu)
  target_include_directories(phx_kernel PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_compile_features(phx_kernel PRIVATE cuda_std_17)
  target_compile_options(phx_kernel PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wno-unknown-pragmas>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
  )
  set_target_properties(phx_kernel PROPERTIES CUDA_SEPARABLE_COMPILATION OFF POSITION_INDEPENDENT_CODE ON)

  target_sources(PHOENIX PRIVATE $<TARGET_OBJECTS:phx_kernel>)
  target_link_libraries(PHOENIX PRIVATE CUDA::cudart)
endif()
