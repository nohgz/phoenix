cmake_minimum_required(VERSION 3.20)

project(PHOENIX LANGUAGES CXX)

# --- Colorized Output ---
string(ASCII 27 ESCAPE)
string(CONCAT PHOENIX_HEAD "[${ESCAPE}[1m${ESCAPE}[38;5;208mPHOENIX${ESCAPE}[0m]")

option(PHX_FORCE_CPU "${PHOENIX_HEAD} Forcing CPU-only compilation, disabling GPU auto-detection" OFF)

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake/Modules"
)
find_package(HDF5 REQUIRED)
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
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets")

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
  add_compile_definitions(PHX_GPU_ENABLED=1)
else()
  if(NOT PHX_FORCE_CPU)
    message(STATUS "${PHOENIX_HEAD} No physical GPU detected. Falling back to CPU acceleration.")
  endif()
  add_compile_definitions(PHX_GPU_ENABLED=0)
endif()

set(HOST_SRC
  src/main.cpp
)

set(CUDA_SRC
  src/device/kernel.cu
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

# Attach device structures exclusively if
# active hardware layer compiled successfully
if(GPU_FOUND)
  target_sources(PHOENIX PRIVATE ${CUDA_SRC})
  target_link_libraries(PHOENIX PRIVATE CUDA::toolkit)

  target_compile_features(PHOENIX PRIVATE
    cuda_std_17
  )

  target_compile_options(PHOENIX PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
  )

  set_target_properties(PHOENIX PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
  )
endif()
