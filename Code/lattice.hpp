#pragma once
// =============================================================================
// Cubic host lattices: BCC (2 atoms per cubic cell) and FCC (4 atoms per cell)
//
// Sites, in a box of nx x ny x nz cubic cells centred on the origin:
//   o + a0 (i, j, k) + a0 f,   o = -a0 (nx, ny, nz) / 2,
//   f in {(0,0,0), (1/2,1/2,1/2)}                              (BCC)
//   f in {(0,0,0), (1/2,1/2,0), (1/2,0,1/2), (0,1/2,1/2)}      (FCC)
// Both contain the cube corners o + a0 Z^3, which snapToMidPlane uses as the
// origin of the atomic planes.
// =============================================================================

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

enum class Lattice { BCC, FCC };

inline const char* latticeName(Lattice l) { return l == Lattice::FCC ? "FCC" : "BCC"; }

// "BCC" / "FCC" (case-sensitive, as the other YAML keywords); false if unknown
inline bool parseLattice(const std::string& s, Lattice& out)
{
    if (s == "BCC") { out = Lattice::BCC; return true; }
    if (s == "FCC") { out = Lattice::FCC; return true; }
    return false;
}

inline int atomsPerCell(Lattice l) { return l == Lattice::FCC ? 4 : 2; }

// Perfect lattice, positions centred at the origin. Atom order: cells (i, j, k),
// k fastest, then the basis sites in the order of the table above.
inline std::vector<std::vector<double>> generateLattice(Lattice l, int nx, int ny, int nz, double a0)
{
    static const double bcc[2][3] = { {0, 0, 0}, {0.5, 0.5, 0.5} };
    static const double fcc[4][3] = { {0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5} };
    const double (*f)[3] = (l == Lattice::FCC) ? fcc : bcc;
    const int nb = atomsPerCell(l);

    std::vector<std::vector<double>> atoms;
    atoms.reserve((size_t)nb * nx * ny * nz);
    const double hx = nx * a0 / 2., hy = ny * a0 / 2., hz = nz * a0 / 2.;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                for (int b = 0; b < nb; b++)
                    atoms.push_back({ i*a0 + f[b][0]*a0 - hx, j*a0 + f[b][1]*a0 - hy, k*a0 + f[b][2]*a0 - hz });
    return atoms;
}

// Spacing of the atomic planes (hkl), h, k, l coprime integers, |n| = sqrt(h^2+k^2+l^2):
//   BCC: a0/|n| if h+k+l is even, else a0/(2|n|)
//   FCC: a0/|n| if h, k, l are all odd, else a0/(2|n|)
inline double latticePlaneSpacing(Lattice l, const long h[3], double a0)
{
    const double norm = std::sqrt(double(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]));
    const bool full = (l == Lattice::FCC)
        ? ((labs(h[0]) % 2) && (labs(h[1]) % 2) && (labs(h[2]) % 2))
        : (labs(h[0] + h[1] + h[2]) % 2 == 0);
    return full ? a0 / norm : a0 / (2. * norm);
}

// Is b = (a0/2) m a translation of the lattice? Only decided for integer m
// (true otherwise): BCC needs m_x, m_y, m_z of the same parity (1/2<111>,
// <100>, ...), FCC needs m_x + m_y + m_z even (1/2<110>, <100>, ...).
inline bool isLatticeTranslation(Lattice l, const std::vector<double>& m)
{
    long h[3];
    for (int i = 0; i < 3; i++) {
        h[i] = std::lround(m[i]);
        if (std::fabs(m[i] - h[i]) > 1e-9) return true;
    }
    if (l == Lattice::FCC) return labs(h[0] + h[1] + h[2]) % 2 == 0;
    return (labs(h[0]) % 2 == labs(h[1]) % 2) && (labs(h[1]) % 2 == labs(h[2]) % 2);
}

// Standard atomic mass (Dalton) of common BCC / FCC metals, 0 if not listed
inline double standardAtomicMass(const std::string& element)
{
    static const struct { const char* name; double mass; } table[] = {
        // BCC
        {"Ta", 180.9479}, {"W", 183.84}, {"Mo", 95.95}, {"Nb", 92.90637}, {"V", 50.9415},
        {"Cr", 51.9961}, {"Fe", 55.845}, {"Li", 6.94}, {"Na", 22.98977}, {"K", 39.0983},
        // FCC
        {"Cu", 63.546}, {"Al", 26.9815385}, {"Ni", 58.6934}, {"Ag", 107.8682}, {"Au", 196.966569},
        {"Pt", 195.084}, {"Pd", 106.42}, {"Pb", 207.2}, {"Ir", 192.217}, {"Rh", 102.9055},
    };
    for (const auto& e : table)
        if (element == e.name) return e.mass;
    return 0.;
}
