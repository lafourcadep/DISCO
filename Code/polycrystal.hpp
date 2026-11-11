#pragma once
// =============================================================================
// Periodic Voronoi polycrystal
//
// Grain g has a seed s_g (lab frame, box [-L/2, L/2)) and a rotation R_g (rows:
// crystal directions of the lab axes, lab = R_g crystal, as in Orientation).
// Its lattice has a site on the seed: lab = s_g + R_g c, c a site of the cubic
// BCC / FCC lattice. An atom belongs to grain g when its nearest seed image is
// s_g itself (periodic Voronoi cell of g).
//
// Filling. A periodic Voronoi cell lies within s_g +/- L/2 on every axis (it is
// closer to s_g than to s_g's own images), so grain g is filled from its lattice
// sites in that box, unwrapped around s_g, and each kept site is wrapped back
// into the box. A grain that crosses a box face is thus one continuous lattice.
// Ties (a site equidistant from two seed images) go to the smaller key
// (grain, image shift), so cells are half-open and no site is kept twice: one
// grain with a commensurate orientation reproduces generateOrientedLattice.
//
// Grain boundaries. Sites of neighbouring grains can come closer than in the
// perfect crystal. Every pair closer than r_min (GB_min_distance x NN distance)
// loses one atom: the one of the higher grain index (or the higher atom index
// in the same grain, where a grain meets its own periodic image). Perfect
// crystal has no such pair, so only GB atoms are affected. The GBs are not
// relaxed: the sample must be relaxed (minimisation, annealing) before use.
// =============================================================================

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "lattice.hpp"

struct Grain {
    double seed[3] = { 0., 0., 0. };   // lab frame, in [-L/2, L/2)
    Orientation o;                     // o.R: rows = crystal directions of the lab axes
    long   atoms = 0;                  // atoms of the grain after GB overlap removal
};

struct PolycrystalSpec {
    bool   active = false;
    double box[3] = { 0., 0., 0. };    // Lx, Ly, Lz (Angstrom)
    int    nGrains = 0;                // random grains (seeds and orientations)
    uint64_t seed = 1;                 // random-number seed
    double gbMinDistance = 0.8;        // fraction of the nearest-neighbour distance
    // explicit grains (Grain_list): centre (lab, box-centred) and orientation;
    // hasCenter / hasOrientation false -> random
    std::vector<Grain> list;
    std::vector<char>  hasCenter, hasOrientation;
};

struct PolycrystalResult {
    std::vector<double> x, y, z;       // atom positions (lab, box [-L/2, L/2))
    std::vector<int>    grain;         // grain index of each atom
    std::vector<Grain>  grains;
    long candidates = 0;               // lattice sites tested
    long removedGB  = 0;               // atoms removed by the GB overlap rule
};

// Deterministic, platform-independent random numbers: the mt19937_64 sequence
// is fixed by the C++ standard, but std::uniform_real_distribution is not, so
// the mapping to [0, 1) is done here (53 random bits).
class PortableRandom {
public:
    explicit PortableRandom(uint64_t seed) : eng_(seed) {}
    double uniform() { return (double)(eng_() >> 11) * (1.0 / 9007199254740992.0); }
private:
    std::mt19937_64 eng_;
};

// Uniformly distributed random rotation (Shoemake), as an Orientation (R only)
Orientation randomOrientation(PortableRandom& rng);

// Bunge Euler angles (degrees) <-> Orientation. R = g^T, g the Bunge matrix
// (crystal coordinates of a sample vector: c = g s).
Orientation orientationFromEuler(double phi1, double Phi, double phi2);
void eulerFromOrientation(const Orientation& o, double& phi1, double& Phi, double& phi2);

// Nearest-neighbour distance of the lattice (FCC a0/sqrt2, BCC a0 sqrt3/2)
double nearestNeighbourDistance(Lattice l, double a0);

// Seeds and orientations from the spec (random ones drawn in grain order)
std::vector<Grain> makeGrains(const PolycrystalSpec& spec);

// Fill the grains with their lattices and apply the GB overlap rule
PolycrystalResult buildPolycrystal(Lattice l, double a0, const PolycrystalSpec& spec,
                                   const std::vector<Grain>& grains);

// Index of the grain whose periodic Voronoi cell holds the lab point x
// (minimum image over all seeds; ties to the smaller grain index)
int grainOf(const std::vector<Grain>& grains, const double box[3], const double x[3]);
