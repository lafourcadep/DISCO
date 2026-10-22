#pragma once
// =============================================================================
// Kroupa pair term — shared by the CUDA kernels (backend_cuda.cu) and the
// OpenMP backend (backend_cpu.cpp), so both evaluate exactly the same formula.
//
// Ref: F. Kroupa, Dislocation dipoles and dislocation loops,
//      J. Phys. Colloques 27 (1966) C3-154.
//
// Displacement at an atom from one loopon (infinitesimal loop of area dS, unit
// glide direction m, unit normal n, Burgers components b), rho = atom - loopon:
//   u_i = -k1 dS b_i / rho^2 [ (1-2nu)/rho (n_i m.rho + m_i n.rho - rho_i m.n)
//                              + 3 rho_i (m.rho)(n.rho) / rho^3 ]
// with k1 = 1 / (8 pi (1-nu)), nu = 0.35.
// =============================================================================

#include <cmath>

#ifdef __CUDACC__
#define KROUPA_HD __host__ __device__ __forceinline__
#else
#define KROUPA_HD inline
#endif

constexpr double KROUPA_NU = 0.35;

// Physical parameters of one dislocation loop (loop-constant)
struct LoopPhys {
    double bx, by, bz;      // Burgers vector components
    double mx, my, mz;      // normalised glide direction m
    double nx, ny, nz;      // normalised loop normal n
    double dS;              // loopon surface element
};

KROUPA_HD double kroupaK1()
{
    return 1.0 / (8.0 * M_PI * (1.0 - KROUPA_NU));
}

// Adds the contribution of one loopon image at separation rho = atom - loopon.
// An atom sitting on the loopon (rho < 1e-10) contributes nothing; the test is
// written as a select rather than a branch so the CPU loop vectorises (the
// arithmetic on accepted pairs is unchanged).
KROUPA_HD void kroupaPair(double rho_x, double rho_y, double rho_z,
                          const LoopPhys& L, double k1, double b_n,
                          double& dux, double& duy, double& duz)
{
    const double rho2 = rho_x*rho_x + rho_y*rho_y + rho_z*rho_z;
    const double rho0 = sqrt(rho2);
    const bool   keep = !(rho0 < 1e-10);        // skip near-zero (atom on loopon)
    const double rho  = keep ? rho0 : 1.0;
    const double r2   = keep ? rho2 : 1.0;
    const double rho3 = r2 * rho;

    const double b_rho = L.mx*rho_x + L.my*rho_y + L.mz*rho_z;   // m . rho
    const double n_rho = L.nx*rho_x + L.ny*rho_y + L.nz*rho_z;   // n . rho
    const double k2    = (1.0 - 2.0*KROUPA_NU) / rho;
    const double fac   = keep ? -k1 * L.dS / r2 : 0.0;
    const double tail  = 3.0 * b_rho * n_rho / rho3;

    dux += fac * L.bx * (k2*(L.nx*b_rho + L.mx*n_rho - rho_x*b_n) + tail*rho_x);
    duy += fac * L.by * (k2*(L.ny*b_rho + L.my*n_rho - rho_y*b_n) + tail*rho_y);
    duz += fac * L.bz * (k2*(L.nz*b_rho + L.mz*n_rho - rho_z*b_n) + tail*rho_z);
}
