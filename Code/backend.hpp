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

// Loopon positions of one loop (SoA) and its physical parameters
struct LoopData {
    std::vector<double> x, y, z;
    LoopPhys phys;
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
#ifdef WITH_CUDA
bool gpuAvailable();                 // a CUDA device can be used
void printGpuInfo();                 // "--- GPU: <name> ---" banner
void kroupaGPU(const AtomArrays& atoms, const std::vector<LoopData>& loops,
               double dx, double dy, double dz, BackendTiming& t);
#endif
