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

static vector<double> bccX, bccY, bccZ;
static void buildBCC(int N)
{
    const double h = N * A0 / 2.;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++)
        for (double f : { 0.0, 0.5 }) {
            bccX.push_back(i*A0 + f*A0 - h);
            bccY.push_back(j*A0 + f*A0 - h);
            bccZ.push_back(k*A0 + f*A0 - h);
        }
}

static void checkSlab(loop& lp, double d, int N, MatterType expect, const char* label);

static void checkLoop(vector<double> m, vector<double> n, double radius, int N, MatterType expect)
{
    double d;
    const vector<double> c = snapToMidPlane({ 1.1, -0.7, 0.4 }, n, A0, N, N, N, d);
    loop lp(0.4, radius, c, m, n, A0, { 1., 1., 1. });
    checkSlab(lp, d, N, expect, "circle ");
}

// Hexagon (regular, optional rounded corners) with the same slab checks
static void checkHexagon(vector<double> m, vector<double> n, double R, double rc,
                         vector<double> vdir, int N, MatterType expect)
{
    double d;
    const vector<double> c = snapToMidPlane({ 1.1, -0.7, 0.4 }, n, A0, N, N, N, d);
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
    const SlabResult sl = prismaticSlab(lp, d, bccX, bccY, bccZ, L, L, L);
    assert(sl.type == expect);
    assert(sl.planes == (int)lround(fabs(bn) / d));
    const double area = lp.getNumberOfLoopons() * lp.getSurface();
    const double nominal = 2. / (A0*A0*A0) * area * fabs(bn);
    const double ns = (double)sl.atoms.size();
    assert(fabs(ns - nominal) < 0.15 * nominal);
    const double tLo = -floor(sl.planes / 2.), tHi = ceil(sl.planes / 2.);
    for (size_t k = 0; k < sl.atoms.size(); k++) {
        const int a = sl.atoms[k];
        const double s = (bccX[a]-c[0])*nh[0] + (bccY[a]-c[1])*nh[1] + (bccZ[a]-c[2])*nh[2];
        assert(s / d >= tLo - 1e-9 && s / d < tHi);
        assert((s > 0.) == (sl.side[k] > 0));
    }
    cout << "  " << label << "m = (" << m0[0] << " " << m0[1] << " " << m0[2] << "), n = ("
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

int main()
{
    const int N = 24;
    buildBCC(N);
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
    checkPolygonGeometry();
    // Glide loop: no matter
    loop g(0.4, 10., { 0., 0., 0. }, { 1, 1, 1 }, { 1, -1, 0 }, A0, { 1., 1., 1. });
    assert(matterType(g) == MatterType::NONE);
    cout << "test_prismatic: OK\n";
    return 0;
}
