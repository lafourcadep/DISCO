#pragma once
// =============================================================================
// Coherent twin lamellae (FCC {111}, BCC {112})
//
// A twin is the mirror image of its parent lattice across a twin plane through
// atoms (both lattices are centrosymmetric, so the mirror is also a 180 deg
// rotation about the plane normal: R_twin = R_parent * halfTurn(n_crystal)).
// A lamella lies between two mirror planes, layers 0 and m along n (spacing d):
// atoms of layers 1 .. m-1 are moved to the twin lattice, layers 0 and m are
// shared. Two mirrors compose to the translation 2 m d n, which must be a
// lattice translation for the matrix on both sides to be one crystal:
// 2 m = 0 mod (period / d) along n, i.e. m a multiple of 3 for FCC {111} and
// BCC {112}. The requested thickness is snapped to the nearest such m.
//
// Construction: in every lamella layer j the twin net is the parent net shifted
// in the plane by one vector t_j (the in-plane lattice is the same), found from
// a reference site: t_j = Mirror(z) - x_ref, x_ref a parent site of layer j and
// z the parent site closest to Mirror(x_ref) (layer -j). Every atom of layer j
// moves by t_j (|t_j| <= a0/sqrt(6) in FCC): one shift per layer keeps the
// mapping one-to-one. No atom is added or removed.
//
// Extent: a lamella spans its grain (polycrystal: bounded by the GBs) or the
// periodic box (single crystal: n must be along a box axis).
// =============================================================================

#include <string>
#include <vector>

#include "lattice.hpp"

struct TwinSpec {                       // as read from the input
    int grain = -1;                     // polycrystal grain (0-based); -1: single crystal
    std::vector<double> n;              // twin plane, crystal frame of the parent
    std::vector<double> center;         // a point on the first twin boundary (box frame)
    double thickness = 0.;              // Angstrom (snapped)
};

struct Twin {
    int    grain = -1;
    long   h[3] = { 0, 0, 0 };          // reduced plane indices (crystal frame)
    double nLab[3] = { 0., 0., 0. };    // unit normal, lab frame
    double O[3] = { 0., 0., 0. };       // parent lattice site on the first TB (lab; unwrapped
                                        //   around the grain seed in a polycrystal)
    double d = 0.;                      // plane spacing
    int    m = 0;                       // layers 0 and m are the two TBs
    int    axis = -1;                   // single crystal: box axis along n (period along n = box)
    int    layersPerBox = 0;            // single crystal: box length / d
    double vLab[3] = { 0., 0., 0. };    // lattice vector with v.n = d (lab)
    std::vector<std::vector<double>> shift;   // t_j, j = 0 .. m (t_0 = t_m = 0)
    Orientation parent, twin;           // crystal frames (twin = parent * half turn)
    long   atomsMoved = 0;
};

// Build a lamella. origin: a parent lattice site (lab: box corner, or grain seed);
// seed: grain seed for unwrapping (nullptr for a single crystal); box: lengths.
// Stops with a message on a non-twin plane or an impossible single-crystal setup.
Twin makeTwin(const TwinSpec& spec, Lattice l, double a0, const Orientation& parent,
              const double origin[3], const double* seed, const double box[3], double& snapShift);

// Layer index of the lab point x (unwrapped as the twin's frame requires):
// 0 .. m = lamella and its TBs; for a single crystal, reduced modulo the box.
long twinLayer(const Twin& t, const double x[3], const double* seed, const double box[3]);

// Move the atoms of the lamella (only those of t.grain in a polycrystal) to the
// twin lattice; atomTwin[i] = index + 1 for moved atoms (0 stays matrix).
void applyTwin(Twin& t, int index, std::vector<std::vector<double>>& atoms,
               const std::vector<int>& atomGrain, std::vector<int>& atomTwin,
               const double* seed, const double box[3]);

// Region of a lab point for twin-aware loops: 0 = matrix, t + 1 = inside lamella t
// (strictly between its TBs), among the twins of the given grain (-1: single crystal).
int twinRegion(const std::vector<Twin>& twins, int grain, const double x[3],
               const double* seed, const double box[3]);
