// =============================================================================
// Kroupa dislocation-loop displacement field — CUDA backend
//
// Parallelism strategy
// --------------------
//   • One GPU thread handles one atom (embarrassingly parallel outer loop).
//   • All 256 threads in a block cooperatively load the current loop's loopon
//     positions into shared memory before computing, eliminating 27×nLp
//     repeated global-memory reads per atom (one read per thread instead of
//     27 × nLp reads per atom from global memory).
//   • One kernel launch per dislocation loop; displacements are accumulated
//     in persistent device arrays d_displX/Y/Z (initialised to zero once).
//   • A tiled fallback kernel is provided for the rare case where the shared
//     memory requirement (3 × nLoopons × 8 bytes) exceeds the GPU limit.
//   • --fp32: single-precision variants of both kernels (see SINGLE PRECISION).
//   • Chunking: a launch can cover only a range of atoms, the host waiting after
//     each one, so that no single GPU wait is long (mvcc on Apple GPUs aborts a
//     wait longer than 8 s). Each atom's sum is unchanged, so are the results.
// The pair term itself is kroupaPairT (kroupa.hpp), shared with the CPU backend.
// =============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <climits>
#include <string>
#include <cuda_runtime.h>

#include "backend.hpp"

using namespace std;

// -----------------------------------------------------------------------------
// CUDA error-checking macro
// -----------------------------------------------------------------------------
#define CUDA_CHECK(call)                                                         \
  do {                                                                           \
    cudaError_t _e = (call);                                                     \
    if (_e != cudaSuccess) {                                                     \
      cerr << "CUDA error at " << __FILE__ << ":" << __LINE__                   \
           << "  " << cudaGetErrorString(_e) << "\n";                            \
      exit(EXIT_FAILURE);                                                        \
    }                                                                            \
  } while (0)

// =============================================================================
// PRIMARY KERNEL  (shared memory — fast path)
//
// Each block cooperatively loads the nLoopons loopon positions into shared
// memory (3 × nLoopons doubles, allocated dynamically at launch).  Every
// thread then iterates over the 27 periodic images and all loopons entirely
// from shared memory, accumulating the displacement in registers.
// =============================================================================
__global__ void kroupa_kernel(
    // Atom positions (Structure of Arrays — coalesced reads)
    const double* __restrict__ d_atomX,
    const double* __restrict__ d_atomY,
    const double* __restrict__ d_atomZ,
    // Accumulated displacements (SoA, initialised to zero before first loop)
    double* d_displX,
    double* d_displY,
    double* d_displZ,
    int nbAtoms,
    // Loopon positions for the current dislocation loop (SoA, in global mem)
    const double* __restrict__ d_lpX,
    const double* __restrict__ d_lpY,
    const double* __restrict__ d_lpZ,
    int nLoopons,
    // Physical parameters of the loop
    LoopPhys L,
    // Supercell dimensions (for 27-image periodic boundary sum)
    double dx, double dy, double dz)
{
    // --- Load loopon positions into shared memory (cooperative) ---
    // Layout: [lpX_0..lpX_{nLp-1} | lpY_0..lpY_{nLp-1} | lpZ_0..lpZ_{nLp-1}]
    extern __shared__ double shLp[];
    double* shX = shLp;
    double* shY = shLp + nLoopons;
    double* shZ = shLp + 2 * nLoopons;

    for (int i = threadIdx.x; i < nLoopons; i += blockDim.x) {
        shX[i] = d_lpX[i];
        shY[i] = d_lpY[i];
        shZ[i] = d_lpZ[i];
    }
    __syncthreads();

    // --- Thread → atom mapping ---
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nbAtoms) return;

    // --- Loop constants: k1 and b.n hoisted outside all loops ---
    const double k1  = kroupaK1();
    const double b_n = L.mx*L.nx + L.my*L.ny + L.mz*L.nz;

    // --- Atom position (registers) ---
    const double ax = d_atomX[tid];
    const double ay = d_atomY[tid];
    const double az = d_atomZ[tid];

    // --- Accumulate displacement in registers ---
    double dux = 0., duy = 0., duz = 0.;

    for (int ii = -1; ii <= 1; ii++) {
    for (int jj = -1; jj <= 1; jj++) {
    for (int kk = -1; kk <= 1; kk++) {

        const double delta_x = dx * ii;
        const double delta_y = dy * jj;
        const double delta_z = dz * kk;

        for (int lp = 0; lp < nLoopons; lp++)
            kroupaPair(ax - (shX[lp] + delta_x), ay - (shY[lp] + delta_y), az - (shZ[lp] + delta_z),
                       L, k1, b_n, dux, duy, duz);
    }}}

    // --- Write accumulated displacement (no race: each thread owns its atom) ---
    d_displX[tid] += dux;
    d_displY[tid] += duy;
    d_displZ[tid] += duz;
}

// =============================================================================
// TILED FALLBACK KERNEL  (static shared memory — for very large loops)
//
// Used when 3 × nLoopons × 8 bytes exceeds the GPU's shared memory per block.
// The host launches this kernel ceil(nLoopons / TILE_SIZE) times, each time
// processing a contiguous tile of loopons.  Displacements are accumulated over
// tiles.
// =============================================================================
static const int TILE_SIZE = 1365;   // 1365 × 3 × 8 = 32,760 bytes < 48 KB

__global__ void kroupa_kernel_tiled(
    const double* __restrict__ d_atomX,
    const double* __restrict__ d_atomY,
    const double* __restrict__ d_atomZ,
    double* d_displX, double* d_displY, double* d_displZ,
    int nbAtoms,
    const double* __restrict__ d_lpX,
    const double* __restrict__ d_lpY,
    const double* __restrict__ d_lpZ,
    int tileOffset, int tileLen,      // process loopons [tileOffset, tileOffset+tileLen)
    LoopPhys L,
    double dx, double dy, double dz)
{
    __shared__ double shX[TILE_SIZE];
    __shared__ double shY[TILE_SIZE];
    __shared__ double shZ[TILE_SIZE];

    for (int i = threadIdx.x; i < tileLen; i += blockDim.x) {
        shX[i] = d_lpX[tileOffset + i];
        shY[i] = d_lpY[tileOffset + i];
        shZ[i] = d_lpZ[tileOffset + i];
    }
    __syncthreads();

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nbAtoms) return;

    const double k1  = kroupaK1();
    const double b_n = L.mx*L.nx + L.my*L.ny + L.mz*L.nz;

    const double ax = d_atomX[tid];
    const double ay = d_atomY[tid];
    const double az = d_atomZ[tid];

    double dux = 0., duy = 0., duz = 0.;

    for (int ii = -1; ii <= 1; ii++) {
    for (int jj = -1; jj <= 1; jj++) {
    for (int kk = -1; kk <= 1; kk++) {
        const double delta_x = dx * ii;
        const double delta_y = dy * jj;
        const double delta_z = dz * kk;

        for (int lp = 0; lp < tileLen; lp++)
            kroupaPair(ax - (shX[lp] + delta_x), ay - (shY[lp] + delta_y), az - (shZ[lp] + delta_z),
                       L, k1, b_n, dux, duy, duz);
    }}}

    d_displX[tid] += dux;
    d_displY[tid] += duy;
    d_displZ[tid] += duz;
}

// =============================================================================
// SINGLE PRECISION  (--fp32)
//
// Apple GPUs have no fp64 hardware (mvcc emulates it), and consumer NVIDIA GPUs
// run fp64 at 1/32-1/64 of fp32. These kernels evaluate the pair term in float:
//   • loopons are uploaded as float offsets from a reference point c of the loop
//     (the loopon centroid), and each atom position is reduced to a - c in double
//     once per thread, then rounded to float: coordinates near the loop, where the
//     field varies fastest, keep ~1e-6 relative precision;
//   • each periodic image's sum over loopons is accumulated in float, and the 27
//     image sums in double.
// =============================================================================
static const int TILE_SIZE_F = 2730;   // 2730 × 3 × 4 = 32,760 bytes

__device__ __forceinline__ void kroupaAtomF(
    double ax_d, double ay_d, double az_d, double cx, double cy, double cz,
    const float* shX, const float* shY, const float* shZ, int nLp,
    const LoopPhysF& L, float k1, float b_n, float dx, float dy, float dz,
    double& ux, double& uy, double& uz)
{
    const float ax = (float)(ax_d - cx);
    const float ay = (float)(ay_d - cy);
    const float az = (float)(az_d - cz);
    for (int ii = -1; ii <= 1; ii++) {
    for (int jj = -1; jj <= 1; jj++) {
    for (int kk = -1; kk <= 1; kk++) {
        const float delta_x = dx * (float)ii;
        const float delta_y = dy * (float)jj;
        const float delta_z = dz * (float)kk;
        float px = 0.f, py = 0.f, pz = 0.f;
        for (int lp = 0; lp < nLp; lp++)
            kroupaPairT<float>(ax - (shX[lp] + delta_x), ay - (shY[lp] + delta_y),
                               az - (shZ[lp] + delta_z), L, k1, b_n, px, py, pz);
        ux += (double)px;
        uy += (double)py;
        uz += (double)pz;
    }}}
}

__global__ void kroupa_kernel_f(
    const double* __restrict__ d_atomX, const double* __restrict__ d_atomY,
    const double* __restrict__ d_atomZ,
    double* d_displX, double* d_displY, double* d_displZ, int nbAtoms,
    const float* __restrict__ d_lpX, const float* __restrict__ d_lpY,
    const float* __restrict__ d_lpZ, int nLoopons,
    LoopPhysF L, float k1, float b_n, double cx, double cy, double cz,
    float dx, float dy, float dz)
{
    extern __shared__ float shLpF[];
    float* shX = shLpF;
    float* shY = shLpF + nLoopons;
    float* shZ = shLpF + 2 * nLoopons;
    for (int i = threadIdx.x; i < nLoopons; i += blockDim.x) {
        shX[i] = d_lpX[i];
        shY[i] = d_lpY[i];
        shZ[i] = d_lpZ[i];
    }
    __syncthreads();

    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nbAtoms) return;
    double ux = 0., uy = 0., uz = 0.;
    kroupaAtomF(d_atomX[tid], d_atomY[tid], d_atomZ[tid], cx, cy, cz,
                shX, shY, shZ, nLoopons, L, k1, b_n, dx, dy, dz, ux, uy, uz);
    d_displX[tid] += ux;
    d_displY[tid] += uy;
    d_displZ[tid] += uz;
}

__global__ void kroupa_kernel_tiled_f(
    const double* __restrict__ d_atomX, const double* __restrict__ d_atomY,
    const double* __restrict__ d_atomZ,
    double* d_displX, double* d_displY, double* d_displZ, int nbAtoms,
    const float* __restrict__ d_lpX, const float* __restrict__ d_lpY,
    const float* __restrict__ d_lpZ, int tileOffset, int tileLen,
    LoopPhysF L, float k1, float b_n, double cx, double cy, double cz,
    float dx, float dy, float dz)
{
    __shared__ float shX[TILE_SIZE_F];
    __shared__ float shY[TILE_SIZE_F];
    __shared__ float shZ[TILE_SIZE_F];
    for (int i = threadIdx.x; i < tileLen; i += blockDim.x) {
        shX[i] = d_lpX[tileOffset + i];
        shY[i] = d_lpY[tileOffset + i];
        shZ[i] = d_lpZ[tileOffset + i];
    }
    __syncthreads();

    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nbAtoms) return;
    double ux = 0., uy = 0., uz = 0.;
    kroupaAtomF(d_atomX[tid], d_atomY[tid], d_atomZ[tid], cx, cy, cz,
                shX, shY, shZ, tileLen, L, k1, b_n, dx, dy, dz, ux, uy, uz);
    d_displX[tid] += ux;
    d_displY[tid] += uy;
    d_displZ[tid] += uz;
}

// =============================================================================
// HOST HELPERS
// =============================================================================
static const int THREADS_PER_BLOCK = 256;

// Device buffers shared by every launch
struct DevBuffers {
    const double *ax, *ay, *az;
    double *ux, *uy, *uz;
    int nbAtoms;
};

// Calls launch(offset, count) over atom ranges of at most `chunkAtoms` atoms
// (all atoms at once when chunkAtoms <= 0), waiting for the GPU after each range
// when chunking.
template <typename F>
static void forAtomChunks(int nbAtoms, int chunkAtoms, F launch)
{
    if (chunkAtoms <= 0 || chunkAtoms >= nbAtoms) {
        launch(0, nbAtoms);
        CUDA_CHECK(cudaGetLastError());
        if (chunkAtoms > 0) CUDA_CHECK(cudaDeviceSynchronize());
        return;
    }
    for (int off = 0; off < nbAtoms; off += chunkAtoms) {
        launch(off, min(chunkAtoms, nbAtoms - off));
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
}

// Atoms per launch for a pair budget and nLp loopons per launch (0 = no chunking)
static int chunkAtomsFor(double chunkPairs, int nLp)
{
    if (chunkPairs <= 0.) return 0;
    long c = (long)(chunkPairs / max(nLp, 1));
    c = max<long>(THREADS_PER_BLOCK, (c / THREADS_PER_BLOCK) * THREADS_PER_BLOCK);
    return (int)min<long>(c, INT_MAX);
}

// One loop, double precision: fast path (all loopons in shared memory) or tiles.
// Returns the number of launches.
static int launchKroupa(const DevBuffers& B,
                        const double* d_lpX, const double* d_lpY, const double* d_lpZ, int nLp,
                        const LoopPhys& L, double dx, double dy, double dz,
                        size_t maxShmem, double chunkPairs)
{
    const size_t shmem = (size_t)3 * nLp * sizeof(double);
    int launches = 0;
    if (shmem <= maxShmem) {
        forAtomChunks(B.nbAtoms, chunkAtomsFor(chunkPairs, nLp), [&](int off, int cnt) {
            kroupa_kernel<<<(cnt + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK, shmem>>>(
                B.ax + off, B.ay + off, B.az + off, B.ux + off, B.uy + off, B.uz + off, cnt,
                d_lpX, d_lpY, d_lpZ, nLp, L, dx, dy, dz);
            launches++;
        });
    } else {
        for (int offset = 0; offset < nLp; offset += TILE_SIZE) {
            const int tileLen = min(TILE_SIZE, nLp - offset);
            forAtomChunks(B.nbAtoms, chunkAtomsFor(chunkPairs, tileLen), [&](int off, int cnt) {
                kroupa_kernel_tiled<<<(cnt + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
                    B.ax + off, B.ay + off, B.az + off, B.ux + off, B.uy + off, B.uz + off, cnt,
                    d_lpX, d_lpY, d_lpZ, offset, tileLen, L, dx, dy, dz);
                launches++;
            });
        }
    }
    return launches;
}

// One loop, single precision (loopons already float offsets from c)
static int launchKroupaF(const DevBuffers& B,
                         const float* d_lpX, const float* d_lpY, const float* d_lpZ, int nLp,
                         const LoopPhysF& L, float k1, float b_n,
                         double cx, double cy, double cz, float dx, float dy, float dz,
                         size_t maxShmem, double chunkPairs)
{
    const size_t shmem = (size_t)3 * nLp * sizeof(float);
    int launches = 0;
    if (shmem <= maxShmem) {
        forAtomChunks(B.nbAtoms, chunkAtomsFor(chunkPairs, nLp), [&](int off, int cnt) {
            kroupa_kernel_f<<<(cnt + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK, shmem>>>(
                B.ax + off, B.ay + off, B.az + off, B.ux + off, B.uy + off, B.uz + off, cnt,
                d_lpX, d_lpY, d_lpZ, nLp, L, k1, b_n, cx, cy, cz, dx, dy, dz);
            launches++;
        });
    } else {
        for (int offset = 0; offset < nLp; offset += TILE_SIZE_F) {
            const int tileLen = min(TILE_SIZE_F, nLp - offset);
            forAtomChunks(B.nbAtoms, chunkAtomsFor(chunkPairs, tileLen), [&](int off, int cnt) {
                kroupa_kernel_tiled_f<<<(cnt + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
                    B.ax + off, B.ay + off, B.az + off, B.ux + off, B.uy + off, B.uz + off, cnt,
                    d_lpX, d_lpY, d_lpZ, offset, tileLen, L, k1, b_n, cx, cy, cz, dx, dy, dz);
                launches++;
            });
        }
    }
    return launches;
}

// =============================================================================
// Device queries
// =============================================================================
bool gpuAvailable()
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

void printGpuInfo()
{
    int devId;
    CUDA_CHECK(cudaGetDevice(&devId));
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, devId));
    cout << "\n--- GPU: " << prop.name << " ---\n";
    cout << "  SMs              : " << prop.multiProcessorCount << "\n";
    cout << "  Global memory    : " << prop.totalGlobalMem / (1 << 20) << " MB\n";
    cout << "  Shared mem/block : " << prop.sharedMemPerBlock / 1024 << " KB\n\n";
}

// =============================================================================
// kroupaGPU — upload atoms, kernel launches per loop, download displacements
// =============================================================================
// Automatic pair budgets per launch on Apple GPUs (mvcc), ~1-2 s per launch:
// measured ~5 M pairs/s in fp64 (emulated) on an M4 Pro.
static const double APPLE_CHUNK_PAIRS_FP64 = 1.0e7;
static const double APPLE_CHUNK_PAIRS_FP32 = 2.0e8;

void kroupaGPU(const AtomArrays& atoms, const vector<LoopData>& loops,
               double dx, double dy, double dz, const GpuOptions& opt, BackendTiming& t)
{
    using Clock = chrono::steady_clock;
    auto seconds = [](Clock::time_point a, Clock::time_point b) {
        return chrono::duration<double>(b - a).count();
    };
    const int nbAtoms = atoms.n;
    const int nbLoops = (int)loops.size();

    int devId;
    CUDA_CHECK(cudaGetDevice(&devId));
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, devId));
    const size_t maxShmem = prop.sharedMemPerBlock;
    const bool appleGPU = string(prop.name).rfind("Apple", 0) == 0;

    double chunkPairs = opt.chunkPairs;
    if (chunkPairs < 0.)
        chunkPairs = appleGPU ? (opt.fp32 ? APPLE_CHUNK_PAIRS_FP32 : APPLE_CHUNK_PAIRS_FP64) : 0.;
    cout << "Precision: " << (opt.fp32 ? "single (fp32 pair term, fp64 sums)" : "double") << "\n";
    if (chunkPairs > 0.)
        cout << "Chunking : <= " << chunkPairs << " atom-loopon pairs per launch, host waits after each"
             << (opt.chunkPairs < 0. ? "  (automatic: Apple GPU)" : "") << "\n";

    // ---- Allocate device memory ----
    Clock::time_point t0 = Clock::now();
    int maxLp = 0;
    for (auto& l : loops) maxLp = max(maxLp, (int)l.x.size());

    double *d_atomX, *d_atomY, *d_atomZ;
    double *d_displX, *d_displY, *d_displZ;
    double *d_lpX = nullptr, *d_lpY = nullptr, *d_lpZ = nullptr;
    float  *d_lpXf = nullptr, *d_lpYf = nullptr, *d_lpZf = nullptr;

    CUDA_CHECK(cudaMalloc(&d_atomX,  nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_atomY,  nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_atomZ,  nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_displX, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_displY, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_displZ, nbAtoms * sizeof(double)));
    if (opt.fp32) {
        CUDA_CHECK(cudaMalloc(&d_lpXf, maxLp * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_lpYf, maxLp * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_lpZf, maxLp * sizeof(float)));
    } else {
        CUDA_CHECK(cudaMalloc(&d_lpX, maxLp * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_lpY, maxLp * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_lpZ, maxLp * sizeof(double)));
    }

    // ---- Upload atom positions ----
    CUDA_CHECK(cudaMemcpy(d_atomX, atoms.x, nbAtoms*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_atomY, atoms.y, nbAtoms*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_atomZ, atoms.z, nbAtoms*sizeof(double), cudaMemcpyHostToDevice));

    // ---- Initialise displacement arrays to zero ----
    CUDA_CHECK(cudaMemset(d_displX, 0, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_displY, 0, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_displZ, 0, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaDeviceSynchronize());
    Clock::time_point t1 = Clock::now();
    t.upload = seconds(t0, t1);

    const DevBuffers B = { d_atomX, d_atomY, d_atomZ, d_displX, d_displY, d_displZ, nbAtoms };

    // ---- Kernel launches per dislocation loop, timed on the host (wait per loop) ----
    cout << "Running GPU kernels...\n";
    t.loopSec.assign(nbLoops, 0.);
    for (int l = 0; l < nbLoops; l++) {
        const Clock::time_point tl = Clock::now();
        const LoopData& lp = loops[l];
        const int nLp = (int)lp.x.size();
        const size_t shmem = (size_t)3 * nLp * (opt.fp32 ? sizeof(float) : sizeof(double));

        cout << "  Loop " << l+1 << "/" << nbLoops
             << ":  " << nLp << " loopons  |  "
             << shmem / 1024.0 << " KB shared mem  "
             << (shmem <= maxShmem ? "[fast path]" : "[tiled]");

        int launches;
        if (opt.fp32) {
            // Reference point c = loopon centroid; loopons as float offsets from c
            double cx = 0., cy = 0., cz = 0.;
            for (int p = 0; p < nLp; p++) { cx += lp.x[p]; cy += lp.y[p]; cz += lp.z[p]; }
            if (nLp > 0) { cx /= nLp; cy /= nLp; cz /= nLp; }
            vector<float> fx(nLp), fy(nLp), fz(nLp);
            for (int p = 0; p < nLp; p++) {
                fx[p] = (float)(lp.x[p] - cx);
                fy[p] = (float)(lp.y[p] - cy);
                fz[p] = (float)(lp.z[p] - cz);
            }
            CUDA_CHECK(cudaMemcpy(d_lpXf, fx.data(), nLp*sizeof(float), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lpYf, fy.data(), nLp*sizeof(float), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lpZf, fz.data(), nLp*sizeof(float), cudaMemcpyHostToDevice));
            const LoopPhys& P = lp.phys;
            const LoopPhysF LF = { (float)P.bx, (float)P.by, (float)P.bz,
                                   (float)P.mx, (float)P.my, (float)P.mz,
                                   (float)P.nx, (float)P.ny, (float)P.nz, (float)P.dS };
            const float k1  = (float)kroupaK1();
            const float b_n = (float)(P.mx*P.nx + P.my*P.ny + P.mz*P.nz);
            launches = launchKroupaF(B, d_lpXf, d_lpYf, d_lpZf, nLp, LF, k1, b_n, cx, cy, cz,
                                     (float)dx, (float)dy, (float)dz, maxShmem, chunkPairs);
        } else {
            CUDA_CHECK(cudaMemcpy(d_lpX, lp.x.data(), nLp*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lpY, lp.y.data(), nLp*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lpZ, lp.z.data(), nLp*sizeof(double), cudaMemcpyHostToDevice));
            launches = launchKroupa(B, d_lpX, d_lpY, d_lpZ, nLp, lp.phys, dx, dy, dz,
                                    maxShmem, chunkPairs);
        }
        CUDA_CHECK(cudaDeviceSynchronize());
        t.loopSec[l] = seconds(tl, Clock::now());
        cout << "  " << launches << " launch" << (launches > 1 ? "es" : "") << ", "
             << t.loopSec[l] << " s\n";
    }
    Clock::time_point t2 = Clock::now();
    t.kernels = seconds(t1, t2);
    cout << "\nGPU kernel time: " << t.kernels * 1e3 << " ms\n\n";

    // ---- Download displacement results ----
    CUDA_CHECK(cudaMemcpy(atoms.ux, d_displX, nbAtoms*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(atoms.uy, d_displY, nbAtoms*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(atoms.uz, d_displZ, nbAtoms*sizeof(double), cudaMemcpyDeviceToHost));

    // ---- Free device memory ----
    cudaFree(d_atomX);  cudaFree(d_atomY);  cudaFree(d_atomZ);
    cudaFree(d_displX); cudaFree(d_displY); cudaFree(d_displZ);
    if (d_lpX)  { cudaFree(d_lpX);  cudaFree(d_lpY);  cudaFree(d_lpZ); }
    if (d_lpXf) { cudaFree(d_lpXf); cudaFree(d_lpYf); cudaFree(d_lpZf); }
    t.download = seconds(t2, Clock::now());
}
