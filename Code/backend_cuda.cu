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
// The pair term itself is kroupaPair (kroupa.hpp), shared with the CPU backend.
// =============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
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
// HOST HELPER — launch the appropriate kernel for one dislocation loop
// =============================================================================
static void launchKroupa(
    // device arrays
    const double* d_atomX, const double* d_atomY, const double* d_atomZ,
    double* d_displX, double* d_displY, double* d_displZ,
    int nbAtoms,
    // device loopon arrays and count
    const double* d_lpX, const double* d_lpY, const double* d_lpZ, int nLp,
    // loop physics
    const LoopPhys& L,
    // box
    double dx, double dy, double dz,
    // GPU properties
    size_t maxShmem,
    // grid/block
    dim3 grid, dim3 block)
{
    size_t shmem = (size_t)3 * nLp * sizeof(double);

    if (shmem <= maxShmem) {
        // Fast path: entire loopon array fits in shared memory
        kroupa_kernel<<<grid, block, shmem>>>(
            d_atomX, d_atomY, d_atomZ,
            d_displX, d_displY, d_displZ, nbAtoms,
            d_lpX, d_lpY, d_lpZ, nLp,
            L, dx, dy, dz);
        CUDA_CHECK(cudaGetLastError());
    } else {
        // Tiled fallback: process loopons in TILE_SIZE-sized chunks
        cout << "    [tiled fallback: " << nLp << " loopons, tile=" << TILE_SIZE << "]\n";
        for (int offset = 0; offset < nLp; offset += TILE_SIZE) {
            int tileLen = min(TILE_SIZE, nLp - offset);
            kroupa_kernel_tiled<<<grid, block>>>(
                d_atomX, d_atomY, d_atomZ,
                d_displX, d_displY, d_displZ, nbAtoms,
                d_lpX, d_lpY, d_lpZ,
                offset, tileLen,
                L, dx, dy, dz);
            CUDA_CHECK(cudaGetLastError());
        }
    }
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
// kroupaGPU — upload atoms, one kernel launch per loop, download displacements
// =============================================================================
void kroupaGPU(const AtomArrays& atoms, const vector<LoopData>& loops,
               double dx, double dy, double dz, BackendTiming& t)
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

    // ---- Allocate device memory ----
    Clock::time_point t0 = Clock::now();
    int maxLp = 0;
    for (auto& l : loops) maxLp = max(maxLp, (int)l.x.size());

    double *d_atomX, *d_atomY, *d_atomZ;
    double *d_displX, *d_displY, *d_displZ;
    double *d_lpX, *d_lpY, *d_lpZ;

    CUDA_CHECK(cudaMalloc(&d_atomX,  nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_atomY,  nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_atomZ,  nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_displX, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_displY, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_displZ, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lpX, maxLp * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lpY, maxLp * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lpZ, maxLp * sizeof(double)));

    // ---- Upload atom positions ----
    CUDA_CHECK(cudaMemcpy(d_atomX, atoms.x, nbAtoms*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_atomY, atoms.y, nbAtoms*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_atomZ, atoms.z, nbAtoms*sizeof(double), cudaMemcpyHostToDevice));

    // ---- Initialise displacement arrays to zero ----
    CUDA_CHECK(cudaMemset(d_displX, 0, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_displY, 0, nbAtoms * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_displZ, 0, nbAtoms * sizeof(double)));
    Clock::time_point t1 = Clock::now();
    t.upload = seconds(t0, t1);

    // ---- Kernel launch configuration ----
    const int threadsPerBlock = 256;
    dim3 block(threadsPerBlock);
    dim3 grid((nbAtoms + threadsPerBlock - 1) / threadsPerBlock);

    // ---- Time the GPU work: ev[l] -> ev[l+1] brackets loop l (loopon upload + kernel) ----
    vector<cudaEvent_t> ev(nbLoops + 1);
    for (auto& e : ev) CUDA_CHECK(cudaEventCreate(&e));
    CUDA_CHECK(cudaEventRecord(ev[0]));

    // ---- One kernel launch per dislocation loop ----
    cout << "Running GPU kernels...\n";
    for (int l = 0; l < nbLoops; l++) {
        const LoopData& lp = loops[l];
        int nLp = (int)lp.x.size();
        size_t shmem = (size_t)3 * nLp * sizeof(double);

        cout << "  Loop " << l+1 << "/" << nbLoops
             << ":  " << nLp << " loopons  |  "
             << shmem / 1024.0 << " KB shared mem  ";
        cout << (shmem <= maxShmem ? "[fast path]\n" : "[tiled]\n");

        // Upload loopon positions
        CUDA_CHECK(cudaMemcpy(d_lpX, lp.x.data(), nLp*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_lpY, lp.y.data(), nLp*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_lpZ, lp.z.data(), nLp*sizeof(double), cudaMemcpyHostToDevice));

        launchKroupa(d_atomX, d_atomY, d_atomZ,
                     d_displX, d_displY, d_displZ, nbAtoms,
                     d_lpX, d_lpY, d_lpZ, nLp,
                     lp.phys, dx, dy, dz,
                     maxShmem, grid, block);
        CUDA_CHECK(cudaEventRecord(ev[l + 1]));
    }

    CUDA_CHECK(cudaEventSynchronize(ev[nbLoops]));
    t.loopSec.assign(nbLoops, 0.);
    for (int l = 0; l < nbLoops; l++) {
        float ms;
        CUDA_CHECK(cudaEventElapsedTime(&ms, ev[l], ev[l + 1]));
        t.loopSec[l] = ms / 1000.;
    }
    float gpuMs;
    CUDA_CHECK(cudaEventElapsedTime(&gpuMs, ev[0], ev[nbLoops]));
    cout << "\nGPU kernel time: " << gpuMs << " ms\n\n";
    for (auto& e : ev) CUDA_CHECK(cudaEventDestroy(e));
    Clock::time_point t2 = Clock::now();
    t.kernels = seconds(t1, t2);

    // ---- Download displacement results ----
    CUDA_CHECK(cudaMemcpy(atoms.ux, d_displX, nbAtoms*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(atoms.uy, d_displY, nbAtoms*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(atoms.uz, d_displZ, nbAtoms*sizeof(double), cudaMemcpyDeviceToHost));

    // ---- Free device memory ----
    cudaFree(d_atomX);  cudaFree(d_atomY);  cudaFree(d_atomZ);
    cudaFree(d_displX); cudaFree(d_displY); cudaFree(d_displZ);
    cudaFree(d_lpX);    cudaFree(d_lpY);    cudaFree(d_lpZ);
    t.download = seconds(t2, Clock::now());
}
