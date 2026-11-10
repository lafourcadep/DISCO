// =============================================================================
// Matter insertion / removal for loops with b.n != 0 — see prismatic.hpp
// =============================================================================

#include "prismatic.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

MatterType matterType(const loop& lp)
{
    const vector<double> m = lp.getM(), n = lp.getN();
    const double mn = m[0]*n[0] + m[1]*n[1] + m[2]*n[2];
    if (fabs(mn) < 1e-6) return MatterType::NONE;
    return mn > 0. ? MatterType::VACANCY : MatterType::INTERSTITIAL;
}

const char* matterName(MatterType t)
{
    switch (t) {
        case MatterType::VACANCY:      return "vacancy";
        case MatterType::INTERSTITIAL: return "interstitial";
        default:                       return "glide";
    }
}

// Minimum image of one component in a periodic box of length L
static inline double minImage(double r, double L)
{
    return r - L * floor(r / L + 0.5);
}

SlabResult prismaticSlab(const loop& lp, double spacing,
                         const vector<double>& x, const vector<double>& y,
                         const vector<double>& z, double lx, double ly, double lz)
{
    SlabResult res;
    res.type = matterType(lp);
    if (res.type == MatterType::NONE) return res;

    const vector<double> m = lp.getM(), n = lp.getN(), c = lp.getCenter();
    const double bmag = lp.getBurgers()[0];          // |b| (compress = 1)
    for (int i = 0; i < 3; i++) res.b[i] = bmag * m[i];
    res.height  = fabs(res.b[0]*n[0] + res.b[1]*n[1] + res.b[2]*n[2]);
    res.spacing = spacing;

    const double ratio = res.height / spacing;
    res.planes = (int)lround(ratio);
    if (res.planes < 1 || fabs(ratio - res.planes) > 1e-6) {
        cerr << "ERROR: |b.n| = " << res.height << " Ang is not a whole number of atomic-plane"
             << " spacings (d = " << spacing << " Ang): this b cannot close on this habit plane\n";
        abort();
    }
    // Slab in units of d, bounded by mid-planes (atoms sit at half-integer s/d)
    const double tLo = -floor(res.planes / 2.);
    const double tHi =  ceil(res.planes / 2.);

    const int nAtoms = (int)x.size();
    for (int a = 0; a < nAtoms; a++) {
        const double rel[3] = { minImage(x[a] - c[0], lx),
                                minImage(y[a] - c[1], ly),
                                minImage(z[a] - c[2], lz) };
        const double s = rel[0]*n[0] + rel[1]*n[1] + rel[2]*n[2];
        const double t = s / spacing;
        if (t < tLo || t >= tHi) continue;
        if (!lp.contains(rel)) continue;
        res.atoms.push_back(a);
        res.side.push_back(s > 0. ? 1 : -1);
    }
    return res;
}
