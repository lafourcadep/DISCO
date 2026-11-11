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

#include <algorithm>
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

// -----------------------------------------------------------------------------
// Oriented crystals. R has as rows the unit crystal directions of the box axes
// X, Y, Z, so a crystal vector c is (R c) in the lab frame. The box spans whole
// periods of the lattice along X, Y and Z (latticePeriod), so it is periodic.
// -----------------------------------------------------------------------------
struct Orientation {
    bool   active = false;
    long   dir[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };   // integer X, Y, Z (crystal frame)
    double R[3][3]   = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
};

inline std::vector<double> rotateToLab(const Orientation& o, const std::vector<double>& c)
{
    if (!o.active || c.empty()) return c;
    std::vector<double> r(3, 0.);
    for (int i = 0; i < 3; i++) r[i] = o.R[i][0]*c[0] + o.R[i][1]*c[1] + o.R[i][2]*c[2];
    return r;
}

inline std::vector<double> rotateToCrystal(const Orientation& o, const std::vector<double>& v)
{
    if (!o.active || v.empty()) return v;
    std::vector<double> c(3, 0.);
    for (int i = 0; i < 3; i++) c[i] = o.R[0][i]*v[0] + o.R[1][i]*v[1] + o.R[2][i]*v[2];
    return c;
}

// Shortest lattice translation along the integer direction d (reduced to coprime
// [uvw]): (a0/2)|uvw| if (a0/2)[uvw] is a translation (FCC: u+v+w even; BCC: u, v,
// w all odd), else a0|uvw|.
inline double latticePeriod(Lattice l, const long d[3], double a0)
{
    long g = 0;
    for (int i = 0; i < 3; i++) { long a = labs(d[i]); while (a) { long t = g % a; g = a; a = t; } }
    if (g == 0) return 0.;
    const std::vector<double> u = { double(d[0] / g), double(d[1] / g), double(d[2] / g) };
    const double norm = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    // isLatticeTranslation is declared below; inline the parity rule here
    const long h0 = labs(d[0] / g), h1 = labs(d[1] / g), h2 = labs(d[2] / g);
    const bool half = (l == Lattice::FCC) ? ((h0 + h1 + h2) % 2 == 0)
                                          : (h0 % 2 == 1 && h1 % 2 == 1 && h2 % 2 == 1);
    return half ? 0.5 * a0 * norm : a0 * norm;
}

// Perfect oriented lattice in the box [-L/2, L/2) (lab frame), with a lattice
// site on the box corner -L/2. Cubic cells covering the rotated box are swept;
// a site is kept if its lab coordinate t (from the corner) is in [-tol, L - tol)
// on every axis, so a site and its periodic image are never both kept.
inline std::vector<std::vector<double>> generateOrientedLattice(Lattice l, const Orientation& o,
                                                               const double L[3], double a0)
{
    static const double bcc[2][3] = { {0, 0, 0}, {0.5, 0.5, 0.5} };
    static const double fcc[4][3] = { {0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5} };
    const double (*f)[3] = (l == Lattice::FCC) ? fcc : bcc;
    const int nb = atomsPerCell(l);
    const double tol = 1e-6;

    // crystal-frame bounding box of the 8 lab corners of [0, L]^3
    double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
    for (int c = 0; c < 8; c++) {
        const double x[3] = { (c & 1) ? L[0] : 0., (c & 2) ? L[1] : 0., (c & 4) ? L[2] : 0. };
        for (int k = 0; k < 3; k++) {
            const double ck = o.R[0][k]*x[0] + o.R[1][k]*x[1] + o.R[2][k]*x[2];
            lo[k] = std::min(lo[k], ck);
            hi[k] = std::max(hi[k], ck);
        }
    }
    long i0[3], i1[3];
    for (int k = 0; k < 3; k++) {
        i0[k] = (long)std::floor(lo[k] / a0) - 1;
        i1[k] = (long)std::ceil(hi[k] / a0) + 1;
    }
    std::vector<std::vector<double>> atoms;
    atoms.reserve((size_t)(nb * L[0] * L[1] * L[2] / (a0 * a0 * a0)) + 16);
    for (long i = i0[0]; i <= i1[0]; i++)
        for (long j = i0[1]; j <= i1[1]; j++)
            for (long k = i0[2]; k <= i1[2]; k++)
                for (int b = 0; b < nb; b++) {
                    const double c[3] = { (i + f[b][0]) * a0, (j + f[b][1]) * a0, (k + f[b][2]) * a0 };
                    double t[3];
                    bool in = true;
                    for (int d = 0; d < 3 && in; d++) {
                        t[d] = o.R[d][0]*c[0] + o.R[d][1]*c[1] + o.R[d][2]*c[2];
                        in = (t[d] >= -tol && t[d] < L[d] - tol);
                    }
                    if (in) atoms.push_back({ t[0] - 0.5 * L[0], t[1] - 0.5 * L[1], t[2] - 0.5 * L[2] });
                }
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

// Is b = (a0/2) m a translation of the lattice? Every BCC / FCC translation has
// an integer m: BCC needs m_x, m_y, m_z of the same parity (1/2<111>, <100>, ...),
// FCC needs m_x + m_y + m_z even (1/2<110>, <100>, ...). A non-integer m
// (Frank a0/3<111>, Shockley a0/6<112>, ...) is never a translation.
inline bool isLatticeTranslation(Lattice l, const std::vector<double>& m)
{
    long h[3];
    for (int i = 0; i < 3; i++) {
        h[i] = std::lround(m[i]);
        if (std::fabs(m[i] - h[i]) > 1e-9) return false;
    }
    if (l == Lattice::FCC) return labs(h[0] + h[1] + h[2]) % 2 == 0;
    return (labs(h[0]) % 2 == labs(h[1]) % 2) && (labs(h[1]) % 2 == labs(h[2]) % 2);
}

// Kind of Burgers vector b = (a0/2) m of a loop of normal n (Miller indices):
//   PERFECT  : b is a lattice translation (any lattice, any n)
//   FRANK    : FCC, b = a0/3<111> (m = 2/3<111>) with n parallel to b: a faulted,
//              sessile loop (intrinsic fault if vacancy, extrinsic if interstitial)
//   SHOCKLEY : FCC, b = a0/6<112> (m = 1/3<112>) on a {111} plane containing b:
//              a glide loop enclosing an intrinsic stacking fault
//   OTHER    : anything else: the cut surface leaves an arbitrary fault
enum class BurgersClass { PERFECT, FRANK, SHOCKLEY, OTHER };

inline const char* burgersClassName(BurgersClass c)
{
    switch (c) {
        case BurgersClass::PERFECT:  return "perfect";
        case BurgersClass::FRANK:    return "Frank";
        case BurgersClass::SHOCKLEY: return "Shockley";
        default:                     return "non-lattice";
    }
}

inline BurgersClass classifyBurgers(Lattice l, const std::vector<double>& m, const std::vector<double>& n)
{
    if (isLatticeTranslation(l, m)) return BurgersClass::PERFECT;
    if (l != Lattice::FCC) return BurgersClass::OTHER;
    long k[3];                                   // k = 3 m, integer for 1/3 fractions
    for (int i = 0; i < 3; i++) {
        k[i] = std::lround(3. * m[i]);
        if (std::fabs(3. * m[i] - k[i]) > 1e-6) return BurgersClass::OTHER;
    }
    const long a0_ = labs(k[0]), a1 = labs(k[1]), a2 = labs(k[2]);
    if (a0_ == 2 && a1 == 2 && a2 == 2) {        // m = 2/3<111>: n must be parallel to m
        const double c[3] = { m[1]*n[2] - m[2]*n[1], m[2]*n[0] - m[0]*n[2], m[0]*n[1] - m[1]*n[0] };
        const double cn = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
        const double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        return (cn < 1e-6 * nn) ? BurgersClass::FRANK : BurgersClass::OTHER;
    }
    const long lo = std::min(a0_, std::min(a1, a2)), hi = std::max(a0_, std::max(a1, a2));
    if (lo == 1 && hi == 2 && a0_ + a1 + a2 == 4) {  // m = 1/3<112>: n in {111}, m.n = 0
        const bool n111 = std::fabs(n[0]) > 0. && std::fabs(std::fabs(n[0]) - std::fabs(n[1])) < 1e-9
                          && std::fabs(std::fabs(n[0]) - std::fabs(n[2])) < 1e-9;
        const double mn = m[0]*n[0] + m[1]*n[1] + m[2]*n[2];
        return (n111 && std::fabs(mn) < 1e-6) ? BurgersClass::SHOCKLEY : BurgersClass::OTHER;
    }
    return BurgersClass::OTHER;
}

// Does a Shockley loop (classifyBurgers == SHOCKLEY) make a stable intrinsic
// stacking fault? Across the cut the +n side moves by -b relative to the -n side
// (u+ - u- = -b). With layers ...X | Y Z... (X below the cut, ABC stacking),
// shifting Y by -b gives Z (intrinsic fault, ABC|BCA) only if -b equals the
// stacking step X -> Y, i.e. if d nhat - b is a lattice translation; the other
// sense of b puts Y on X (AA-type stacking), which is not a metastable fault.
// The step is the same between any two adjacent {111} layers, so the answer
// does not depend on which mid-plane holds the loop. Reverse m (or n) to switch.
inline bool shockleyFaultIsStable(const std::vector<double>& m, const std::vector<double>& n)
{
    const double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    // in units of a0/2: d nhat = (a0/sqrt3) n/|n| -> (2/sqrt3) n/|n|, b -> m
    const double k = 2. / (std::sqrt(3.) * nn);
    return isLatticeTranslation(Lattice::FCC, { k*n[0] - m[0], k*n[1] - m[1], k*n[2] - m[2] });
}

// Displacement of the copy of an interstitial Frank-loop atom relative to its
// source (FCC). The source layer sits next to the cut; v (= side * b, normal to
// the {111} habit plane, |v| = d111) alone would stack the copy straight over
// its source (AA). Adding the in-plane Shockley vector s (a0/6<112>, s.nhat = 0)
// for which v + s is NOT a lattice translation puts the inserted layer on the
// third stacking position: ...A [C] B C... (extrinsic fault). (With the other
// sense of s, v + s is a translation, i.e. the copy would repeat the layer on
// the far side of the cut.) out = v + s.
inline void frankCopyShift(const double v[3], const double nhat[3], double a0, double out[3])
{
    static const int perm[3][3] = { {0, 1, 2}, {1, 2, 0}, {2, 0, 1} };
    for (int p = 0; p < 3; p++)
        for (int sg = 0; sg < 8; sg++) {
            double s[3];
            const double base[3] = { 1., 1., 2. };
            for (int i = 0; i < 3; i++)
                s[perm[p][i]] = ((sg >> i) & 1 ? -1. : 1.) * base[i] * a0 / 6.;
            if (std::fabs(s[0]*nhat[0] + s[1]*nhat[1] + s[2]*nhat[2]) > 1e-9 * a0) continue;
            const std::vector<double> w = { 2. * (v[0] + s[0]) / a0, 2. * (v[1] + s[1]) / a0,
                                            2. * (v[2] + s[2]) / a0 };
            if (isLatticeTranslation(Lattice::FCC, w)) continue;
            for (int i = 0; i < 3; i++) out[i] = v[i] + s[i];
            return;
        }
    for (int i = 0; i < 3; i++) out[i] = v[i];   // not reached for a {111} normal
}

// Standard atomic mass (Dalton) and room-temperature lattice of common BCC /
// FCC metals; mass 0 if not listed
struct ElementData { const char* name; double mass; Lattice lattice; };

inline const ElementData* findElement(const std::string& element)
{
    static const ElementData table[] = {
        // BCC
        {"Ta", 180.9479, Lattice::BCC}, {"W", 183.84, Lattice::BCC}, {"Mo", 95.95, Lattice::BCC},
        {"Nb", 92.90637, Lattice::BCC}, {"V", 50.9415, Lattice::BCC}, {"Cr", 51.9961, Lattice::BCC},
        {"Fe", 55.845, Lattice::BCC}, {"Li", 6.94, Lattice::BCC}, {"Na", 22.98977, Lattice::BCC},
        {"K", 39.0983, Lattice::BCC}, {"Rb", 85.4678, Lattice::BCC}, {"Cs", 132.90545, Lattice::BCC},
        {"Ba", 137.327, Lattice::BCC}, {"Eu", 151.964, Lattice::BCC},
        // FCC
        {"Cu", 63.546, Lattice::FCC}, {"Al", 26.9815385, Lattice::FCC}, {"Ni", 58.6934, Lattice::FCC},
        {"Ag", 107.8682, Lattice::FCC}, {"Au", 196.966569, Lattice::FCC}, {"Pt", 195.084, Lattice::FCC},
        {"Pd", 106.42, Lattice::FCC}, {"Pb", 207.2, Lattice::FCC}, {"Ir", 192.217, Lattice::FCC},
        {"Rh", 102.9055, Lattice::FCC}, {"Ca", 40.078, Lattice::FCC}, {"Sr", 87.62, Lattice::FCC},
        {"Ce", 140.116, Lattice::FCC}, {"Yb", 173.045, Lattice::FCC}, {"Th", 232.0377, Lattice::FCC},
    };
    for (const auto& e : table)
        if (element == e.name) return &e;
    return nullptr;
}

inline double standardAtomicMass(const std::string& element)
{
    const ElementData* e = findElement(element);
    return e ? e->mass : 0.;
}
