/**
 * @file kernel.cu
 * @author C. Noga (github.com/nohgz)
 * @brief The CUDA ray tracer
 *
 * GNU GPL v3.0
 *
 * Important information about the method.
 *
 * - The grid is scaled so that it fills up the GPU, then each thread consumes pixels from a job counter
 *   until there is no more in the list. This is done to ensure that all warps stay busy throughout the integration.
 * - The job list is a morton z ordered curve to ensure spatial locality of pixels.
 * - The job list is processed in chunks so that a single launch stays short
 *   enough for display watchdogs, and so the host can update a progress bar.
 */

#include <cuda_runtime.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <atomic>
#include <chrono>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "core/numerics/headers.hpp"

#include "device/kernel.cuh"
#include "core/geometry/headers.hpp"
#include "core/model/headers.hpp"
#include "core/render/headers.hpp"

#ifndef GPU_BLOCK
  #define GPU_BLOCK 128  // threads per block
#endif
#ifndef GPU_MINBLOCKS
  #define GPU_MINBLOCKS 3  // __launch_bounds__ min blocks per SM. caps registers at 65536 / (GPU_BLOCK * GPU_MINBLOCKS)
#endif
#ifndef GPU_STACK
  #define GPU_STACK 4096  // default per-thread stack limit (bytes)
#endif
#ifndef GPU_CHUNK
  #define GPU_CHUNK 262144  // pixels per kernel launch
#endif

namespace device {
  __device__ int g_dbg_count;

  /**
   * @brief The active lanes take consecutive numbers from the global counter using one atomic.
   */
  __device__ static inline int grab(int* counter) {
    const unsigned mask = __activemask();
    const int lane = threadIdx.x & 31;
    const int leader = __ffs(mask) - 1;
    const int rank = __popc(mask & ((1u << lane) - 1u));
    int base = 0;
    if (lane == leader) base = atomicAdd(counter, __popc(mask));
    base = __shfl_sync(mask, base, leader);
    return base + rank;
  }

  /**
   * @brief One thread per pixel. Evaluates are from the job list [j0, j1) until it is empty.
   */
  __global__ void __launch_bounds__(GPU_BLOCK, GPU_MINBLOCKS) k_render(
    const __grid_constant__ Model model,
    const int* jobs,
    int j0,
    int j1,
    int* counter,
    FreqOut* res,
    unsigned long long* tot
  ) {
    long long steps = 0, nplasma = 0, bad = 0;
    for (;;) {
      const int n = j0 + grab(counter);
      if (n >= j1) break;
      const int idx = jobs[n];
      const int i = idx / model.uparams.img_height, j = idx % model.uparams.img_height;
      double a0, b0;
      core::camera::pixel_impact(i, j, model.uparams.img_width, model.uparams.img_height, model.uparams.fov_x, model.uparams.fov_y, &a0, &b0);
      const double alpha = a0 * model.cos_pa + b0 * model.sin_pa;
      const double beta = -a0 * model.sin_pa + b0 * model.cos_pa;
      const core::pixel::RunStats ps = core::pixel::render_pixel(&model, alpha, beta, &res[(size_t)idx * model.uparams.nfreqs]);
      if (ps.bad && atomicAdd(&g_dbg_count, 1) < 8)
        printf("Non-finite segment operator: pixel (%d,%d) alpha %g beta %g (%ld segments dropped)\n", i, j, alpha, beta,
               ps.bad);
      steps += ps.steps;
      nplasma += ps.nplasma;
      bad += ps.bad;
    }
    atomicAdd(&tot[0], (unsigned long long)steps);
    atomicAdd(&tot[1], (unsigned long long)nplasma);
    atomicAdd(&tot[2], (unsigned long long)bad);
  }

  /**
   * @brief GPU decode of the GRMHD dump, where we map one thread to one cell.
   */
  __global__ void __launch_bounds__(128)
      k_decode(const __grid_constant__ Model m, const double* raw, int leaf0, int nl, scalar* prim) {
    const size_t t = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (size_t)nl * m.cells) return;
    const int q = (int)(t / m.cells), cell = (int)(t % m.cells), leaf = leaf0 + q;
    core::bhac::convert2prim(&m, m.lb + 3 * leaf, m.dxc + 3 * leaf, raw + (size_t)q * m.nwini * m.cells, cell,
                             prim + ((size_t)leaf * m.cells + cell) * core::bhac::NPRIM);
  }

  // -- host helpers --

  static void check(cudaError_t e, const char* file, int line) {
    if (e != cudaSuccess) {
      throw std::runtime_error(std::string(cudaGetErrorString(e)) + " (" + file + ":" + std::to_string(line) + ")");
    }
  }

#define CK(call) check((call), __FILE__, __LINE__)

  static double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  /**
   * @brief Device allocations that are freed together when the owner goes away.
   */
  class DeviceAllocs {
   public:
    DeviceAllocs() {}
    DeviceAllocs(const DeviceAllocs&) = delete;
    DeviceAllocs& operator=(const DeviceAllocs&) = delete;
    ~DeviceAllocs() {
      for (void* p : ptrs_) cudaFree(p);
    }
    template <class T>
    T* alloc(size_t bytes) {
      void* d = NULL;
      CK(cudaMalloc(&d, bytes ? bytes : 1));
      ptrs_.push_back(d);
      return (T*)d;
    }
    template <class T>
    T* upload(const T* host, size_t bytes) {
      T* d = alloc<T>(bytes);
      if (bytes) CK(cudaMemcpy(d, host, bytes, cudaMemcpyHostToDevice));
      return d;
    }

   private:
    std::vector<void*> ptrs_;
  };

  static std::thread* g_warm = NULL;
  static double g_decode_s = 0.;

  void warmup_start() {
    if (!g_warm) g_warm = new std::thread([] { cudaFree(0); });
  }

  static void warmup_wait() {
    if (g_warm) {
      g_warm->join();
      delete g_warm;
      g_warm = NULL;
    }
  }

  void count_progress(int n, void* user) { *(std::atomic<int>*)user += n; }

  /**
   * @brief Reads a run-time override (positive integer) from the environment.
   */
  static int env_int(const char* name, int fallback) {
    const char* env = getenv(name);
    return (env && atoi(env) > 0) ? atoi(env) : fallback;
  }

  /**
   * @brief Persistent-thread grid: fills the machine once, but no more threads than `limit` jobs.
   */
  static int persistent_grid(int sms, int limit) {
    int bpsm = 0;
    CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&bpsm, k_render, GPU_BLOCK, 0));
    if (bpsm < 1) bpsm = 1;
    int grid = bpsm * sms;
    if (grid * GPU_BLOCK > limit) grid = (limit + GPU_BLOCK - 1) / GPU_BLOCK;  // small images: no idle threads
    return grid < 1 ? 1 : grid;
  }

  // -- model upload --

  /**
   * @brief Streams the dump's leaf blocks to the device through two pinned buffers and two streams (the file read of
   *        chunk k+1 overlaps the copy and decode of chunk k) and decodes them there into d->prim.
   */
  static void decode_to_device(DeviceAllocs& al, const Model* m, const char* dump, Model* d) {
    const auto t0 = std::chrono::steady_clock::now();
    const size_t leaf_doubles = (size_t)m->nwini * m->cells;
    const long stag = (long)(m->nx[0] + 1) * (m->nx[1] + 1) * (m->nx[2] + 1) * m->nws * 8;
    int per = (int)((32u << 20) / (leaf_doubles * sizeof(double)));  // about 32 mb per buffer
    if (per < 1) per = 1;
    if (per > m->nleafs) per = m->nleafs;

    d->prim = al.alloc<scalar>((size_t)m->nleafs * m->cells * core::bhac::NPRIM * sizeof(scalar));
    double *h_raw[2] = {NULL, NULL}, *d_raw[2] = {NULL, NULL};
    cudaStream_t st[2] = {NULL, NULL};
    cudaEvent_t ev[2] = {NULL, NULL};
    FILE* f = NULL;
    auto cleanup = [&]() {
      if (f) fclose(f);
      for (int b = 0; b < 2; b++) {
        if (st[b]) cudaStreamSynchronize(st[b]);
        if (ev[b]) cudaEventDestroy(ev[b]);
        if (st[b]) cudaStreamDestroy(st[b]);
        if (h_raw[b]) cudaFreeHost(h_raw[b]);
        if (d_raw[b]) cudaFree(d_raw[b]);
      }
    };
    try {
      for (int b = 0; b < 2; b++) {
        CK(cudaMallocHost((void**)&h_raw[b], (size_t)per * leaf_doubles * sizeof(double)));
        CK(cudaMalloc((void**)&d_raw[b], (size_t)per * leaf_doubles * sizeof(double)));
        CK(cudaStreamCreate(&st[b]));
        CK(cudaEventCreate(&ev[b]));
        CK(cudaEventRecord(ev[b], st[b]));
      }
      f = fopen(dump, "rb");
      if (!f) {
        throw std::runtime_error(std::string("cannot open ") + dump);
      }
      int nbuf = 0;
      for (int l0 = 0; l0 < m->nleafs; l0 += per, nbuf ^= 1) {
        const int nl = (l0 + per < m->nleafs) ? per : m->nleafs - l0;
        CK(cudaEventSynchronize(ev[nbuf]));
        for (int q = 0; q < nl; q++) {
          if (fread(h_raw[nbuf] + (size_t)q * leaf_doubles, sizeof(double), leaf_doubles, f) != leaf_doubles) {
            throw std::runtime_error("short read in the data, leaf " + std::to_string(l0 + q));
          }
          fseek(f, stag, SEEK_CUR);
        }
        CK(cudaMemcpyAsync(d_raw[nbuf], h_raw[nbuf], (size_t)nl * leaf_doubles * sizeof(double), cudaMemcpyHostToDevice,
                           st[nbuf]));
        CK(cudaEventRecord(ev[nbuf], st[nbuf]));
        const size_t cells = (size_t)nl * m->cells;
        k_decode<<<(unsigned)((cells + 127) / 128), 128, 0, st[nbuf]>>>(*d, d_raw[nbuf], l0, nl, d->prim);
        CK(cudaGetLastError());
      }
      CK(cudaStreamSynchronize(st[0]));
      CK(cudaStreamSynchronize(st[1]));
    } catch (...) {
      cleanup();
      throw;
    }
    cleanup();
    g_decode_s = seconds_since(t0);
  }

  /**
   * @brief Shallow copy of the Model, then every table to the device and the pointers patched (Model is plain data).
   */
  static void upload_model(DeviceAllocs& al, const Model* m, const char* dump, Model* d) {
    const size_t n3 = 3 * (size_t)m->nleafs;
    *d = *m;
    d->lb = al.upload(m->lb, n3 * sizeof(scalar));
    d->dxc = al.upload(m->dxc, n3 * sizeof(scalar));
    if (!m->prim) {
      decode_to_device(al, m, dump, d);
    } else {
      g_decode_s = 0.;
      d->prim = al.upload(m->prim, (size_t)m->nleafs * m->cells * core::bhac::NPRIM * sizeof(scalar));
    }
    for (int L = 0; L < m->uparams.nlevels && L < 16; L++) {
      const size_t n = (size_t)m->level_dims[L][0] * m->level_dims[L][1] * m->level_dims[L][2];
      d->level_map[L] = al.upload(m->level_map[L], n * sizeof(int));
    }
  }

  // -- render --

  void render(
    const Model* model,
    const char* dump,
    const int* jobs,
    int njobs,
    FreqOut* res,
    Info* info,
    Progress progress,
    void* user
  ) {
    const auto t_up = std::chrono::steady_clock::now();
    warmup_wait();
    memset(info, 0, sizeof *info);
    const int chunk = env_int("PHX_CHUNK", GPU_CHUNK), stack = env_int("PHX_STACK", GPU_STACK);
    const size_t nres = (size_t)model->uparams.img_width * model->uparams.img_height * model->uparams.nfreqs;
    unsigned long long tot[3] = {0, 0, 0};

    cudaDeviceProp prop;
    CK(cudaGetDeviceProperties(&prop, 0));
    snprintf(info->device, sizeof info->device, "%s (sm_%d%d)", prop.name, prop.major, prop.minor);
    info->sms = prop.multiProcessorCount;
    CK(cudaDeviceSetLimit(cudaLimitStackSize, (size_t)stack));

    DeviceAllocs al;
    Model d;
    upload_model(al, model, dump, &d);
    int* d_jobs = al.upload(jobs, (size_t)njobs * sizeof(int));
    int* d_counter = al.alloc<int>(sizeof(int));
    FreqOut* d_res = al.upload(res, nres * sizeof(FreqOut));
    unsigned long long* d_tot = al.upload(tot, sizeof tot);
    info->upload_s = seconds_since(t_up);
    info->decode_s = g_decode_s;

    // setup launch geometry & persistent grid
    cudaFuncAttributes fa;
    CK(cudaFuncGetAttributes(&fa, k_render));
    info->regs = fa.numRegs;
    info->local_bytes = (int)fa.localSizeBytes;
    const int grid = persistent_grid(prop.multiProcessorCount, chunk);
    info->blocks = grid;
    info->block = GPU_BLOCK;

    // chunk by chunk raymarching
    cudaEvent_t e0, e1;
    CK(cudaEventCreate(&e0));
    CK(cudaEventCreate(&e1));
    float ms_total = 0.f;
    try {
      for (int j0 = 0; j0 < njobs; j0 += chunk) {
        const int j1 = (j0 + chunk < njobs) ? j0 + chunk : njobs;
        CK(cudaMemset(d_counter, 0, sizeof(int)));
        CK(cudaEventRecord(e0, 0));
        k_render<<<grid, GPU_BLOCK>>>(d, d_jobs, j0, j1, d_counter, d_res, d_tot);
        CK(cudaGetLastError());
        CK(cudaEventRecord(e1, 0));
        CK(cudaEventSynchronize(e1));
        float ms = 0.f;
        CK(cudaEventElapsedTime(&ms, e0, e1));
        ms_total += ms;
        if (progress) progress(j1 - j0, user);
      }
      info->kernel_s = ms_total * 1e-3;

      // results
      CK(cudaMemcpy(res, d_res, nres * sizeof(FreqOut), cudaMemcpyDeviceToHost));
      CK(cudaMemcpy(tot, d_tot, sizeof tot, cudaMemcpyDeviceToHost));
      info->steps = (long)tot[0];
      info->nplasma = (long)tot[1];
      info->bad = (long)tot[2];
    } catch (...) {
      cudaEventDestroy(e0);
      cudaEventDestroy(e1);
      throw;
    }
    cudaEventDestroy(e0);
    cudaEventDestroy(e1);
  }
}
