#pragma once
// =============================================================================
// Matter insertion / removal for loops with b.n != 0 (prismatic or mixed loops)
//
// A Volterra loop whose Burgers vector b has a component normal to its cut
// surface opens a gap (interstitial loop) or makes the two lips overlap
// (vacancy loop) by b.n_hat across the cut. The displacement field alone
// cannot do this: without adding or removing atoms the loop is only a strain
// field, and it relaxes away. This module selects, for each such loop, the
// slab of atoms inside the contour that must be duplicated or deleted, as in
// Atomsk (opt_disloc.f90, "loop").
//
// Sign convention (the one of the Kroupa field in kroupa.hpp, b = |b| m):
//   across the cut, from the -n side to the +n side,  u+ - u- = -b.
//   m.n > 0  -> the lips overlap by |b.n_hat|  -> VACANCY loop (atoms removed)
//   m.n < 0  -> a gap of |b.n_hat| opens       -> INTERSTITIAL loop (atoms added)
//
// Slab: atoms at signed distance s = (x - c).n_hat from the cut (c on a
// mid-plane, see snapToMidPlane) with s/d in [-floor(k/2), ceil(k/2)), d the
// atomic-plane spacing of the family n and k = |b.n_hat| / d (an integer, since
// b is a lattice vector). The slab bounds are mid-planes, so exactly k atomic
// planes are taken whatever the rounding.
//   vacancy      : slab atoms are deleted before the displacement sum.
//   interstitial : every slab atom gets a copy at the same lattice site, whose
//                  displacement is u(x) + sign(s) b: the copy lands on the
//                  other lip of the cut, and originals + copies fill the gap
//                  with perfect crystal (b is a lattice translation).
//                  FCC Frank loop (b = a0/3<111>, not a translation): main.cpp
//                  adds the in-plane Shockley shift of frankCopyShift
//                  (lattice.hpp), so the inserted {111} layer takes the third
//                  stacking position (extrinsic fault) instead of AA stacking.
//                  A vacancy Frank loop needs nothing: removing one layer
//                  leaves an intrinsic fault.
// =============================================================================

#include <vector>

#include "loop.hpp"

enum class MatterType { NONE, VACANCY, INTERSTITIAL };

// Matter type of a loop from b = |b| m and n (NONE for a glide loop, b.n = 0)
MatterType matterType(const loop& lp);

const char* matterName(MatterType t);

struct SlabResult {
    MatterType type = MatterType::NONE;
    int    planes   = 0;       // k: atomic planes in the slab
    double height   = 0.;      // |b.n_hat| (Angstrom)
    double spacing  = 0.;      // d (Angstrom)
    double b[3]     = {0., 0., 0.};   // Burgers vector |b| m (Angstrom)
    std::vector<int> atoms;    // indices of the slab atoms
    std::vector<int> side;     // sign of s for each slab atom (+1 / -1)
};

// Select the slab atoms of loop lp among the atoms (x, y, z) of a periodic
// orthorhombic box (lx, ly, lz) centred on the origin. spacing is the
// atomic-plane spacing d of the loop normal family (from snapToMidPlane).
// Aborts if |b.n_hat| is not a whole number of plane spacings.
SlabResult prismaticSlab(const loop& lp, double spacing,
                         const std::vector<double>& x, const std::vector<double>& y,
                         const std::vector<double>& z, double lx, double ly, double lz);
