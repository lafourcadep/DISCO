#pragma once
// =============================================================================
// Compute backends for the Kroupa displacement sum.
//
// Both backends take the same inputs and accumulate into the same host arrays:
//   for every atom, every dislocation loop, the 27 periodic images of the
//   supercell and every loopon of the loop  (kroupaPair in kroupa.hpp).
//   backend_cpu.cpp  — OpenMP over atoms, SIMD over loopons (always built)
//   backend_cuda.cu  — one CUDA thread per atom (built with -DWITH_CUDA)
// =============================================================================

#include <string>
#include <vector>

#include "kroupa.hpp"

// One loop for the backends: loopon positions (SoA) and its physical parameters
// (Method: loopons), or its contour polyline and fan apex (Method: line,
// barnett.hpp; the Burgers vector is (bx mx, by my, bz mz) from phys)
struct LoopData {
    std::vector<double> x, y, z;
    LoopPhys phys;
    bool line = false;
    std::vector<double> vx, vy, vz;    // polyline vertices, counter-clockwise about n
    double cx = 0., cy = 0., cz = 0.;  // fan apex (loop centre)
};

// Atom positions (SoA, read-only) and accumulated displacements (zeroed by the caller)
struct AtomArrays {
    const double *x, *y, *z;
    double *ux, *uy, *uz;
    int n;
};

// Wall-clock seconds of the backend phases; loopSec[l] = loop l alone
struct BackendTiming {
    double upload = 0., kernels = 0., download = 0.;
    std::vector<double> loopSec;
};

// ---- CPU (OpenMP) ----
int  cpuThreads();
void kroupaCPU(const AtomArrays& atoms, const std::vector<LoopData>& loops,
               double dx, double dy, double dz, BackendTiming& t);

// ---- GPU (CUDA) ----
// fp32       : single-precision kernel (pair term in float, coordinates taken
//              relative to the loop centre, per-image partial sums added in double)
// chunkPairs : atom-loopon pairs (each with its 27 images) per kernel launch; the
//              host waits after each launch. < 0: automatic (on for Apple GPUs run
//              through mvcc, whose runtime aborts a GPU wait longer than 8 s; off
//              otherwise), 0: off (one launch per loop or tile)
struct GpuOptions {
    bool   fp32       = false;
    double chunkPairs = -1.;
};
#ifdef WITH_CUDA
bool gpuAvailable();                 // a CUDA device can be used
void printGpuInfo();                 // "--- GPU: <name> ---" banner
bool gpuLineFp64();                  // the build has the fp64 line kernel (not under mvcc)
void kroupaGPU(const AtomArrays& atoms, const std::vector<LoopData>& loops,
               double dx, double dy, double dz, const GpuOptions& opt, BackendTiming& t);
#endif
