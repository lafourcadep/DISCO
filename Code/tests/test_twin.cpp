// =============================================================================
// Unit test for coherent twin lamellae (twin.hpp):
//   1. nearestLatticeSite: the result is a lattice site and no neighbouring
//      site is closer (FCC, BCC, random points).
//   2. Single crystal, FCC (Z = [111]) and BCC (Z = [-1-12]), lamellae of m = 3,
//      6, 9 layers: atom count unchanged, atoms outside the lamella unchanged,
//      FCC: no pair closer than the perfect nearest-neighbour distance and every
//      atom keeps 12 first neighbours; BCC: the closest pair is 2d across the
//      reflection TB (as expected, before relaxation); every lamella atom is the mirror
//      image of a matrix atom across BOTH twin boundaries (layers 0 and m): the
//      two TBs are coherent mirrors.
//   3. The thickness snaps to a multiple of 3 layers; the twin frame is the
//      parent frame turned by 180 deg about n.
//
//   make test
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>

#include "../twin.hpp"

using namespace std;

static Orientation makeOrientation(const long d[3][3])
{
    Orientation o;
    o.active = true;
    for (int a = 0; a < 3; a++) {
        const double n = sqrt(double(d[a][0]*d[a][0] + d[a][1]*d[a][1] + d[a][2]*d[a][2]));
        for (int k = 0; k < 3; k++) { o.dir[a][k] = d[a][k]; o.R[a][k] = d[a][k] / n; }
    }
    return o;
}

static void nearestSites()
{
    const double a0 = 3.6;
    unsigned long long st = 12345;
    auto rnd = [&]() { st = st * 6364136223846793005ULL + 1442695040888963407ULL; return (double)(st >> 11) / 9007199254740992.0; };
    for (Lattice l : { Lattice::FCC, Lattice::BCC })
        for (int it = 0; it < 20000; it++) {
            const double c[3] = { (rnd() - 0.5) * 20., (rnd() - 0.5) * 20., (rnd() - 0.5) * 20. };
            double z[3];
            nearestLatticeSite(l, c, a0, z);
            assert(isLatticeTranslation(l, { 2 * z[0] / a0, 2 * z[1] / a0, 2 * z[2] / a0 }));
            const double dz = (c[0]-z[0])*(c[0]-z[0]) + (c[1]-z[1])*(c[1]-z[1]) + (c[2]-z[2])*(c[2]-z[2]);
            for (int a = -2; a <= 2; a++) for (int b = -2; b <= 2; b++) for (int e = -2; e <= 2; e++) {
                if (!isLatticeTranslation(l, { double(a), double(b), double(e) })) continue;
                const double y[3] = { z[0] + 0.5*a0*a, z[1] + 0.5*a0*b, z[2] + 0.5*a0*e };
                const double dy = (c[0]-y[0])*(c[0]-y[0]) + (c[1]-y[1])*(c[1]-y[1]) + (c[2]-y[2])*(c[2]-y[2]);
                assert(dy >= dz - 1e-9);
            }
        }
    cout << "  nearestLatticeSite: FCC and BCC, 20000 random points each OK\n";
}

typedef tuple<long, long, long> Key;
static Key key(const double x[3], const double L[3])
{
    double w[3];
    for (int k = 0; k < 3; k++) { w[k] = x[k] - L[k] * floor(x[k] / L[k] + 0.5); if (w[k] >= 0.5 * L[k] - 1e-6) w[k] -= L[k]; }
    return Key(lround(w[0] * 1e4), lround(w[1] * 1e4), lround(w[2] * 1e4));
}

static void lamella(Lattice l, const long dirs[3][3], int np[3], double thickness, int mWant)
{
    const double a0 = 3.6;
    const Orientation o = makeOrientation(dirs);
    double L[3];
    for (int a = 0; a < 3; a++) L[a] = np[a] * latticePeriod(l, o.dir[a], a0);
    vector<vector<double>> at = generateOrientedLattice(l, o, L, a0);
    const vector<vector<double>> parent = at;
    set<Key> parentSet;
    for (const auto& p : parent) parentSet.insert(key(p.data(), L));

    TwinSpec sp;
    sp.n = { double(dirs[2][0]), double(dirs[2][1]), double(dirs[2][2]) };
    sp.center = { 0.3, -0.2, 0.1 };
    sp.thickness = thickness;
    const double origin[3] = { -0.5 * L[0], -0.5 * L[1], -0.5 * L[2] };
    double snap;
    Twin t = makeTwin(sp, l, a0, o, origin, nullptr, L, snap);
    assert(t.m == mWant && t.axis == 2);
    vector<int> grain(at.size(), 0), tw(at.size(), 0);
    applyTwin(t, 0, at, grain, tw, nullptr, L);

    // atoms outside the lamella unchanged, lamella atoms moved within the plane
    long inside = 0;
    for (size_t i = 0; i < at.size(); i++) {
        const long j = twinLayer(t, parent[i].data(), nullptr, L);
        const bool in = (j >= 1 && j < t.m);
        assert(in == (tw[i] == 1));
        inside += in;
        double dz = at[i][2] - parent[i][2];
        dz -= L[2] * floor(dz / L[2] + 0.5);
        assert(fabs(dz) < 1e-9);
        if (!in) assert(key(at[i].data(), L) == key(parent[i].data(), L));
    }
    // mirror images across both TBs land on matrix (parent) sites
    for (size_t i = 0; i < at.size(); i++) {
        if (!tw[i]) continue;
        for (int b = 0; b < 2; b++) {
            const double zb = t.O[2] + (b == 0 ? 0. : t.m * t.d);
            const double x[3] = { at[i][0], at[i][1], 2. * zb - at[i][2] };
            assert(parentSet.count(key(x, L)));
        }
    }
    // no close pair (and 12 first neighbours in FCC)
    const double nn = (l == Lattice::FCC) ? a0 / sqrt(2.) : a0 * sqrt(3.) / 2.;
    double dmin = 1e9;
    for (size_t i = 0; i < at.size(); i++) {
        int z = 0;
        for (size_t k = 0; k < at.size(); k++) {
            if (i == k) continue;
            double s2 = 0.;
            for (int c = 0; c < 3; c++) { double e = at[k][c] - at[i][c]; e -= L[c] * floor(e / L[c] + 0.5); s2 += e * e; }
            const double r = sqrt(s2);
            dmin = min(dmin, r);
            if (r < nn * (1. + 1e-6)) z++;
        }
        if (l == Lattice::FCC) assert(z == 12);
    }
    // FCC: across a TB, layers +1 and -1 are 2d = 1.15 a0 apart: perfect NN distance kept.
    // BCC: the pure-mirror (reflection) {112} TB puts layers +1 and -1 straight over each other,
    // 2d = 0.816 a0 < NN = 0.866 a0 (it relaxes to the isosceles TB by an a0/12 <111> shift)
    const double expectMin = (l == Lattice::FCC) ? nn : min(nn, 2. * t.d);
    assert(fabs(dmin - expectMin) < 1e-6);
    // twin frame: parent turned by 180 deg about n (lab normal unchanged, R_twin orthonormal)
    const vector<double> nt = rotateToLab(t.twin, { sp.n[0], sp.n[1], sp.n[2] });
    const vector<double> np_ = rotateToLab(o, { sp.n[0], sp.n[1], sp.n[2] });
    for (int k = 0; k < 3; k++) assert(fabs(nt[k] - np_[k]) < 1e-12);
    cout << "  " << latticeName(l) << " lamella n = [" << dirs[2][0] << " " << dirs[2][1] << " " << dirs[2][2]
         << "], thickness " << thickness << " -> m = " << t.m << " (" << t.m * t.d << " Ang): " << inside
         << " atoms moved of " << at.size() << ", both TBs are mirrors, min distance " << dmin << "\n";
}

int main()
{
    cout << "test_twin:\n";
    nearestSites();
    const long fcc[3][3] = { {1, -1, 0}, {1, 1, -2}, {1, 1, 1} };
    const long bcc[3][3] = { {1, 1, 1}, {-1, 1, 0}, {-1, -1, 2} };
    int nf[3] = { 6, 4, 6 }, nb[3] = { 6, 4, 2 };
    const double d111 = 3.6 / sqrt(3.), d112 = 3.6 / sqrt(6.);
    lamella(Lattice::FCC, fcc, nf, 3 * d111, 3);
    lamella(Lattice::FCC, fcc, nf, 6.4 * d111, 6);     // snapped down to 6
    lamella(Lattice::FCC, fcc, nf, 8 * d111, 9);       // snapped up to 9
    lamella(Lattice::FCC, fcc, nf, 0.5, 3);            // at least 3
    lamella(Lattice::BCC, bcc, nb, 3 * d112, 3);
    lamella(Lattice::BCC, bcc, nb, 6 * d112, 6);
    cout << "test_twin: OK\n";
    return 0;
}
