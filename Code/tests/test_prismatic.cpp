// =============================================================================
// Unit test for prismatic loops (prismatic.hpp):
//   1. Sign convention of the Kroupa field: across the cut of a loop with b.n != 0,
//      (u+ - u-).n_hat = -b.n_hat, so m.n > 0 is a vacancy loop and m.n < 0 an
//      interstitial loop, as matterType() assumes.
//   2. The slab holds |b.n|/d atomic planes inside the contour: atom count close to
//      rho * A * |b.n|, all atoms within the slab bounds, sides consistent.
//   3. Mixed loop (1/2[111] on (110)): one plane, same checks.
//   4. Hexagonal loops (sharp and rounded corners): same checks; polygon geometry
//      (square polygon == rectangle, exact contour area, corner-radius clamp).
//   5. FCC lattice: prismatic 1/2[110] on (110) (2 planes), mixed 1/2[101] on (111)
//      (1 plane), a[001] on (001) (2 planes), same checks.
//   6. FCC faulted loops: Frank a0/3<111> loops (1 plane), interstitial copies
//      stacked as an extrinsic fault (frankCopyShift); classifyBurgers and the
//      sense of Shockley loops (shockleyFaultIsStable).
//
//   make test
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../kroupa.hpp"
#include "../loop.hpp"
#include "../prismatic.hpp"

using namespace std;

static const double A0 = 3.3024;

// Normal displacement jump of the loop's own field (no periodic images) across its
// centre, measured at c +/- h n_hat
static double jumpAcrossCut(const loop& lp, double h)
{
    const vector<double> b = lp.getBurgers(), m = lp.getM(), n = lp.getN(), c = lp.getCenter();
    const LoopPhys L = { b[0], b[1], b[2], m[0], m[1], m[2], n[0], n[1], n[2], lp.getSurface() };
    const double k1 = kroupaK1(), b_n = m[0]*n[0] + m[1]*n[1] + m[2]*n[2];
    double un[2];
    for (int side = 0; side < 2; side++) {
        const double sg = side == 0 ? 1. : -1.;
        const double p[3] = { c[0] + sg*h*n[0], c[1] + sg*h*n[1], c[2] + sg*h*n[2] };
        double ux = 0., uy = 0., uz = 0.;
        for (int k = 0; k < lp.getNumberOfLoopons(); k++) {
            const vector<double> q = lp.getOneLoopon(k);
            kroupaPair(p[0]-q[0], p[1]-q[1], p[2]-q[2], L, k1, b_n, ux, uy, uz);
        }
        un[side] = ux*n[0] + uy*n[1] + uz*n[2];
    }
    return un[0] - un[1];
}

// sites of the current lattice (generateLattice), atom density rho
static vector<double> latX, latY, latZ;
static Lattice lat = Lattice::BCC;
static double rho = 0.;
static void buildLattice(Lattice l, int N)
{
    lat = l;
    rho = atomsPerCell(l) / (A0*A0*A0);
    latX.clear(); latY.clear(); latZ.clear();
    for (const auto& p : generateLattice(l, N, N, N, A0)) {
        latX.push_back(p[0]); latY.push_back(p[1]); latZ.push_back(p[2]);
    }
}

static void checkSlab(loop& lp, double d, int N, MatterType expect, const char* label);

static void checkLoop(vector<double> m, vector<double> n, double radius, int N, MatterType expect)
{
    double d;
    const vector<double> c = snapToMidPlane({ 1.1, -0.7, 0.4 }, n, A0, N, N, N, d, lat);
    loop lp(0.4, radius, c, m, n, A0, { 1., 1., 1. });
    checkSlab(lp, d, N, expect, "circle ");
}

// Hexagon (regular, optional rounded corners) with the same slab checks
static void checkHexagon(vector<double> m, vector<double> n, double R, double rc,
                         vector<double> vdir, int N, MatterType expect)
{
    double d;
    const vector<double> c = snapToMidPlane({ 1.1, -0.7, 0.4 }, n, A0, N, N, N, d, lat);
    loop lp(0.4, regularPolygon(6, R, vdir, m, n), rc, c, m, n, A0, { 1., 1., 1. }, vdir);
    assert(lp.getShape() == LoopShape::POLYGON && lp.getNumberOfVertices() == 6);
    checkSlab(lp, d, N, expect, rc > 0. ? "hexagon (rounded) " : "hexagon ");
}

static void checkSlab(loop& lp, double d, int N, MatterType expect, const char* label)
{
    const vector<double> m0 = lp.getM(), n0 = lp.getN(), c = lp.getCenter();
    lp.generateSurface();
    lp.generateLoopons();
    const double bmag = lp.getBurgers()[0];
    const vector<double> mh = lp.getM(), nh = lp.getN();
    const double bn = bmag * (mh[0]*nh[0] + mh[1]*nh[1] + mh[2]*nh[2]);

    // 1. sign of the jump (h small compared with the radius, several loopon steps)
    const double J = jumpAcrossCut(lp, 2.0);
    assert(fabs(J + bn) < 0.15 * fabs(bn));
    assert(matterType(lp) == expect);

    // 2. slab
    const double L = N * A0;
    const SlabResult sl = prismaticSlab(lp, d, latX, latY, latZ, L, L, L);
    assert(sl.type == expect);
    assert(sl.planes == (int)lround(fabs(bn) / d));
    const double area = lp.getNumberOfLoopons() * lp.getSurface();
    const double nominal = rho * area * fabs(bn);
    const double ns = (double)sl.atoms.size();
    assert(fabs(ns - nominal) < 0.15 * nominal);
    const double tLo = -floor(sl.planes / 2.), tHi = ceil(sl.planes / 2.);
    for (size_t k = 0; k < sl.atoms.size(); k++) {
        const int a = sl.atoms[k];
        const double s = (latX[a]-c[0])*nh[0] + (latY[a]-c[1])*nh[1] + (latZ[a]-c[2])*nh[2];
        assert(s / d >= tLo - 1e-9 && s / d < tHi);
        assert((s > 0.) == (sl.side[k] > 0));
    }
    cout << "  " << latticeName(lat) << " " << label << "m = (" << m0[0] << " " << m0[1] << " " << m0[2] << "), n = ("
         << n0[0] << " " << n0[1] << " " << n0[2] << "): " << matterName(sl.type) << ", jump " << J << " Ang (-b.n = "
         << -bn << "), " << sl.planes << " plane(s), " << sl.atoms.size() << " slab atoms (rho A |b.n| = "
         << nominal << ")\n";
}

// Polygon geometry: a square polygon equals the (rounded) rectangle point by point,
// the loopon area tends to the exact contour area, and the corner radius is clamped
// to the apothem of a regular hexagon (which then becomes its inscribed circle).
static void checkPolygonGeometry()
{
    const vector<double> m = { 1, -1, 1 }, n = { 1, 1, 0 }, c = { 0., 0., 0. }, one = { 1., 1., 1. };
    const double h = 10.;
    for (double rc : { 0., 3. }) {
        loop rect(0.4, h, h, rc, c, m, n, A0, one);
        // vertices (+-h, +-h) along m_hat and p_hat = m_hat x n_hat
        const vector<double> mh = rect.getM(), ph = rect.getP();
        vector< vector<double> > sq;
        for (auto sgn : { make_pair(1., 1.), make_pair(-1., 1.), make_pair(-1., -1.), make_pair(1., -1.) }) {
            vector<double> v(3);
            for (int i = 0; i < 3; i++) v[i] = h * (sgn.first * mh[i] + sgn.second * ph[i]);
            sq.push_back(v);
        }
        loop poly(0.4, sq, rc, c, m, n, A0, one);
        unsigned seed = 12345;
        int mismatch = 0;
        for (int k = 0; k < 20000; k++) {
            double rel[3];
            for (int i = 0; i < 3; i++) {
                seed = seed * 1103515245u + 12345u;
                rel[i] = ((seed >> 8) / double(1u << 24) - 0.5) * 3. * h;
            }
            if (rect.contains(rel) != poly.contains(rel)) mismatch++;
        }
        assert(mismatch == 0);
        assert(fabs(rect.getContourArea() - poly.getContourArea()) < 1e-9);
    }

    // Regular hexagon: area 3 sqrt(3)/2 R^2, rounded: area(Q) + perimeter(Q) r + pi r^2
    const double R = 15.;
    for (double rc : { 0., 2., 1e9 }) {
        loop hx(0.1, regularPolygon(6, R, {}, m, n), rc, c, m, n, A0, one);
        hx.generateSurface();
        hx.generateLoopons();
        const double ap = R * sqrt(3.) / 2., r = hx.getCornerRadius();
        const double exact = (rc > 1e8) ? M_PI * ap * ap
                           : 1.5 * sqrt(3.) * R * R - r * r * (6. / sqrt(3.) - M_PI) * 1.0;
        const double grid = hx.getNumberOfLoopons() * hx.getSurface();
        assert(fabs(hx.getContourArea() - exact) < 1e-6 * exact);
        assert(fabs(grid - exact) < 0.01 * exact);
        if (rc > 1e8) assert(fabs(r - ap) < 1e-9);
        cout << "  hexagon R = " << R << ", corner " << r << ": area " << hx.getContourArea()
             << " (loopons " << grid << ")\n";
    }
    // Vertex_dir: the first vertex is along it
    const vector< vector<double> > v = regularPolygon(6, R, { 1, -1, 0 }, { 1, 1, 1 }, { 1, 1, 1 });
    assert(fabs(v[0][0] - R / sqrt(2.)) < 1e-9 && fabs(v[0][1] + R / sqrt(2.)) < 1e-9 && fabs(v[0][2]) < 1e-9);
}

// Is p a site of the FCC lattice of generateLattice (N^3 cells, corner o = -A0 N/2)?
static bool isFccSite(const double p[3], int N)
{
    vector<double> q(3);
    for (int i = 0; i < 3; i++) q[i] = (p[i] + 0.5 * A0 * N) / (0.5 * A0);
    return isLatticeTranslation(Lattice::FCC, { (double)lround(q[0]), (double)lround(q[1]), (double)lround(q[2]) })
        && fabs(q[0] - lround(q[0])) < 1e-6 && fabs(q[1] - lround(q[1])) < 1e-6 && fabs(q[2] - lround(q[2])) < 1e-6;
}

// Interstitial Frank loop: each copy (perfect-lattice frame, before the field) sits
// at the height of the layer X across the cut, so it must be neither on an X site
// (repeat of X) nor straight over its source (AA): the third stacking position.
static void checkFrankCopies(vector<double> m, vector<double> n, int N)
{
    double d;
    const vector<double> c = snapToMidPlane({ 1.1, -0.7, 0.4 }, n, A0, N, N, N, d, lat);
    loop lp(0.4, 15., c, m, n, A0, { 1., 1., 1. });
    lp.generateSurface();
    lp.generateLoopons();
    assert(classifyBurgers(Lattice::FCC, m, n) == BurgersClass::FRANK);
    const double L = N * A0;
    const SlabResult sl = prismaticSlab(lp, d, latX, latY, latZ, L, L, L);
    assert(sl.type == MatterType::INTERSTITIAL && sl.planes == 1 && !sl.atoms.empty());
    const vector<double> nh = lp.getN();
    const double nhat[3] = { nh[0], nh[1], nh[2] };
    for (size_t k = 0; k < sl.atoms.size(); k++) {
        const int a = sl.atoms[k];
        const double v[3] = { sl.side[k]*sl.b[0], sl.side[k]*sl.b[1], sl.side[k]*sl.b[2] };
        double w[3];
        frankCopyShift(v, nhat, A0, w);
        const double sv = (w[0]-v[0])*nhat[0] + (w[1]-v[1])*nhat[1] + (w[2]-v[2])*nhat[2];
        assert(fabs(sv) < 1e-9);                                   // in-plane Shockley shift
        assert(fabs(sqrt((w[0]-v[0])*(w[0]-v[0]) + (w[1]-v[1])*(w[1]-v[1]) + (w[2]-v[2])*(w[2]-v[2]))
                    - A0 / sqrt(6.)) < 1e-9);
        const double p[3]  = { latX[a] + w[0], latY[a] + w[1], latZ[a] + w[2] };
        const double pv[3] = { p[0] - v[0], p[1] - v[1], p[2] - v[2] };   // copy brought to its source layer
        assert(!isFccSite(p, N));      // not a repeat of the layer X across the cut
        assert(!isFccSite(pv, N));     // not straight over its source (AA)
    }
    cout << "  FCC Frank interstitial m = (" << m[0] << " " << m[1] << " " << m[2] << "), n = (" << n[0] << " "
         << n[1] << " " << n[2] << "): " << sl.atoms.size() << " copies on the third stacking position\n";
}

int main()
{
    const int N = 24;
    buildLattice(Lattice::BCC, N);
    cout << "test_prismatic:\n";
    checkLoop({ 1, 1, 1 },   { 1, 1, 1 }, 18., N, MatterType::VACANCY);       // pure prismatic
    checkLoop({ -1, -1, -1 }, { 1, 1, 1 }, 18., N, MatterType::INTERSTITIAL);
    checkLoop({ 1, -1, 1 },  { 1, -1, 0 }, 18., N, MatterType::VACANCY);       // mixed, {110}
    checkLoop({ -1, 1, -1 }, { 1, -1, 0 }, 18., N, MatterType::INTERSTITIAL);
    checkLoop({ 0, 0, 2 },   { 0, 0, 1 }, 18., N, MatterType::VACANCY);        // a[001]: 2 planes
    // Hexagonal prismatic loops on (111), vertices along <110> (sides along <112>)
    checkHexagon({ -1, -1, -1 }, { 1, 1, 1 }, 18., 0., { 1, -1, 0 }, N, MatterType::INTERSTITIAL);
    checkHexagon({ 1, 1, 1 },    { 1, 1, 1 }, 18., 3., { 1, -1, 0 }, N, MatterType::VACANCY);
    checkHexagon({ -1, -1, -1 }, { 1, 1, 1 }, 18., 3., {},           N, MatterType::INTERSTITIAL);
    // FCC: perfect 1/2<110> and a<100> Burgers vectors
    buildLattice(Lattice::FCC, N);
    checkLoop({ 1, 1, 0 },   { 1, 1, 0 }, 18., N, MatterType::VACANCY);       // pure prismatic, 2 planes
    checkLoop({ -1, -1, 0 }, { 1, 1, 0 }, 18., N, MatterType::INTERSTITIAL);
    checkLoop({ 1, 0, 1 },   { 1, 1, 1 }, 18., N, MatterType::VACANCY);       // mixed, {111}: 1 plane
    checkLoop({ -1, 0, -1 }, { 1, 1, 1 }, 18., N, MatterType::INTERSTITIAL);
    checkLoop({ 0, 0, 2 },   { 0, 0, 1 }, 18., N, MatterType::VACANCY);       // a[001]: 2 planes
    checkHexagon({ -1, -1, 0 }, { 1, 1, 0 }, 18., 3., {}, N, MatterType::INTERSTITIAL);
    assert(!isLatticeTranslation(Lattice::FCC, { 1, 1, 1 }) && isLatticeTranslation(Lattice::FCC, { 1, -1, 0 }));
    assert(!isLatticeTranslation(Lattice::BCC, { 1, 1, 0 }) && isLatticeTranslation(Lattice::BCC, { 1, -1, 1 }));
    // FCC faulted loops: Frank a0/3<111> (1 plane of d111 = a0/sqrt3)
    const double t3 = 2. / 3.;
    checkLoop({ t3, t3, t3 },    { 1, 1, 1 }, 18., N, MatterType::VACANCY);
    checkLoop({ -t3, -t3, -t3 }, { 1, 1, 1 }, 18., N, MatterType::INTERSTITIAL);
    checkFrankCopies({ -t3, -t3, -t3 }, { 1, 1, 1 }, N);
    checkFrankCopies({ t3, -t3, t3 },   { -1, 1, -1 }, N);
    checkFrankCopies({ t3, t3, -t3 },   { -1, -1, 1 }, N);
    const double o3 = 1. / 3.;
    assert(!isLatticeTranslation(Lattice::FCC, { t3, t3, t3 }));
    assert(classifyBurgers(Lattice::FCC, { 1, -1, 0 }, { 1, 1, 1 }) == BurgersClass::PERFECT);
    assert(classifyBurgers(Lattice::FCC, { t3, t3, t3 }, { 1, 1, 1 }) == BurgersClass::FRANK);
    assert(classifyBurgers(Lattice::FCC, { t3, t3, t3 }, { 1, -1, 0 }) == BurgersClass::OTHER);   // not on {111}
    assert(classifyBurgers(Lattice::FCC, { o3, o3, -2*o3 }, { 1, 1, 1 }) == BurgersClass::SHOCKLEY);
    assert(classifyBurgers(Lattice::FCC, { o3, o3, -2*o3 }, { 1, -1, 1 }) == BurgersClass::OTHER); // b.n != 0
    assert(classifyBurgers(Lattice::FCC, { 1, 1, 1 }, { 1, 1, 1 }) == BurgersClass::OTHER);
    assert(classifyBurgers(Lattice::BCC, { t3, t3, t3 }, { 1, 1, 1 }) == BurgersClass::OTHER);
    // Shockley sense (checked with OVITO CNA/DXA: the stable sense gives an HCP fault layer pair)
    assert(!shockleyFaultIsStable({ o3, o3, -2*o3 },  { 1, 1, 1 }));
    assert( shockleyFaultIsStable({ -o3, -o3, 2*o3 }, { 1, 1, 1 }));
    assert( shockleyFaultIsStable({ o3, o3, -2*o3 },  { -1, -1, -1 }));
    // Each {111} plane has 3 stable and 3 unstable Shockley vectors
    for (const vector<double>& n : { vector<double>{ 1, 1, 1 }, vector<double>{ -1, 1, 1 }, vector<double>{ 1, -1, 1 } }) {
        int stable = 0, all = 0;
        for (int p = 0; p < 3; p++)
            for (int sg = 0; sg < 8; sg++) {
                vector<double> m(3);
                const double base[3] = { o3, o3, 2*o3 };
                for (int i = 0; i < 3; i++) m[(i + p) % 3] = ((sg >> i) & 1 ? -1. : 1.) * base[i];
                if (classifyBurgers(Lattice::FCC, m, n) != BurgersClass::SHOCKLEY) continue;
                all++;
                stable += shockleyFaultIsStable(m, n);
            }
        assert(all == 6 && stable == 3);
    }
    buildLattice(Lattice::BCC, N);
    checkPolygonGeometry();
    // Glide loop: no matter
    loop g(0.4, 10., { 0., 0., 0. }, { 1, 1, 1 }, { 1, -1, 0 }, A0, { 1., 1., 1. });
    assert(matterType(g) == MatterType::NONE);
    cout << "test_prismatic: OK\n";
    return 0;
}
