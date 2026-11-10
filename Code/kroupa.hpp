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

// Physical parameters of one dislocation loop (loop-constant). R = double for
// the reference kernels, float for the single-precision GPU kernel (--fp32).
template <typename R>
struct LoopPhysT {
    R bx, by, bz;      // Burgers vector components
    R mx, my, mz;      // normalised glide direction m
    R nx, ny, nz;      // normalised loop normal n
    R dS;              // loopon surface element
};
using LoopPhys  = LoopPhysT<double>;
using LoopPhysF = LoopPhysT<float>;

KROUPA_HD double kroupaK1()
{
    return 1.0 / (8.0 * M_PI * (1.0 - KROUPA_NU));
}

KROUPA_HD double kroupaSqrt(double x) { return sqrt(x); }
KROUPA_HD float  kroupaSqrt(float x)  { return sqrtf(x); }

// Adds the contribution of one loopon image at separation rho = atom - loopon.
// An atom sitting on the loopon (rho < 1e-10) contributes nothing; the test is
// written as a select rather than a branch so the CPU loop vectorises (the
// arithmetic on accepted pairs is unchanged).
// Every constant is cast to R, so the float version has no double operation
// (Apple GPUs emulate fp64); for R = double the arithmetic is unchanged.
template <typename R>
KROUPA_HD void kroupaPairT(R rho_x, R rho_y, R rho_z,
                           const LoopPhysT<R>& L, R k1, R b_n,
                           R& dux, R& duy, R& duz)
{
    const R rho2 = rho_x*rho_x + rho_y*rho_y + rho_z*rho_z;
    const R rho0 = kroupaSqrt(rho2);
    const bool keep = !(rho0 < R(1e-10));       // skip near-zero (atom on loopon)
    const R rho  = keep ? rho0 : R(1.0);
    const R r2   = keep ? rho2 : R(1.0);
    const R rho3 = r2 * rho;

    const R b_rho = L.mx*rho_x + L.my*rho_y + L.mz*rho_z;   // m . rho
    const R n_rho = L.nx*rho_x + L.ny*rho_y + L.nz*rho_z;   // n . rho
    const R k2    = R(1.0 - 2.0*KROUPA_NU) / rho;
    const R fac   = keep ? -k1 * L.dS / r2 : R(0.0);
    const R tail  = R(3.0) * b_rho * n_rho / rho3;

    dux += fac * L.bx * (k2*(L.nx*b_rho + L.mx*n_rho - rho_x*b_n) + tail*rho_x);
    duy += fac * L.by * (k2*(L.ny*b_rho + L.my*n_rho - rho_y*b_n) + tail*rho_y);
    duz += fac * L.bz * (k2*(L.nz*b_rho + L.mz*n_rho - rho_z*b_n) + tail*rho_z);
}

KROUPA_HD void kroupaPair(double rho_x, double rho_y, double rho_z,
                          const LoopPhys& L, double k1, double b_n,
                          double& dux, double& duy, double& duz)
{
    kroupaPairT<double>(rho_x, rho_y, rho_z, L, k1, b_n, dux, duy, duz);
}
