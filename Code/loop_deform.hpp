#pragma once
// =============================================================================
// Box deformation from inserted dislocation loops (Kroupa / Orowan approach)
//
// Physics
// -------
//   Plastic distortion:  β[i,j] = Σ_k  b_k[i] * n_k[j] * A_k / V_box
//   Deformation gradient: F = I + β
//   Polar decomposition:  F = R · U   (SVD via eigendecomposition of F^T·F)
//   Apply pure stretch:   H' = U · H,   r' = U · r
//
// Ref: Tadmor & Miller, Modeling Materials (2011), §2.5
//      Cottrell, Dislocations and Plastic Flow in Crystals (1953)
// =============================================================================

#include <array>
#include <vector>

// 3×3 matrix stored as M[row][col]
using Mat3 = std::array<std::array<double,3>,3>;

// --------------------------------------------------------------------------
// Descriptor for one inserted dislocation loop
// --------------------------------------------------------------------------
struct LoopDesc {
    double b[3];   // Burgers vector components [Å]        (direction + magnitude)
    double n[3];   // Glide-plane normal (normalised internally if not unit)
    double A;      // Total swept area [Å²]  =  nLoopons * dS
};

// --------------------------------------------------------------------------
// computeU
//
//   Accumulates the plastic distortion tensor β from all loops, builds
//   F = I + β, and returns the pure right-stretch tensor U obtained via
//   polar decomposition (SVD through eigendecomposition of F^T·F).
//
//   V_box : simulation box volume [Å³]
//
//   Prints a warning if det(U) departs from 1 by more than 0.01 %
//   (expected to be zero for pure shear loops where b·n = 0).
// --------------------------------------------------------------------------
Mat3 computeU(const std::vector<LoopDesc>& loops, double V_box);

// --------------------------------------------------------------------------
// applyLoopDeformation  (convenience wrapper)
//
//   Computes U (see above), then:
//     H'   = U · H        (updates H in-place)
//     r'_k = U · r_k      (updates each atom position in-place)
//
//   Returns U for downstream use (printing, further checks, …).
// --------------------------------------------------------------------------
Mat3 applyLoopDeformation(
    const std::vector<LoopDesc>& loops,
    Mat3&                        H,
    std::vector<double>&         atomX,
    std::vector<double>&         atomY,
    std::vector<double>&         atomZ,
    int                          nAtoms);
