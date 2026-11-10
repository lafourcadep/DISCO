// =============================================================================
// Unit test for the line-integral loop field (barnett.hpp) and loop::polyline:
//   1. Agreement with a direct quadrature of the Kroupa kernel over the cut
//      (7-point Gauss rule on subdivided fan triangles of the same polyline), for
//      glide, prismatic and mixed loops of every shape.
//   2. Jump across the cut: u+ - u- = -b inside the contour, 0 outside.
//   3. Invariance: fan apex moved inside the contour, starting vertex rotated.
//   4. Kroupa loopon sum (loop::generateLoopons) converges to it as da -> 0.
//   5. Far field: a small loop matches one loopon of the same area.
//   6. polyline: counter-clockwise about n, area of the contour, chord sagitta.
//
//   make test
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../barnett.hpp"
#include "../loop.hpp"

using namespace std;
typedef vector<double> V;

static const double A0 = 3.3024;
static const V ONE = { 1., 1., 1. };

static V cross(const V& a, const V& b)
{
    return { a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] };
}
static double dot(const V& a, const V& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static double maxAbs3(const double a[3], const double b[3])
{
    return fmax(fabs(a[0]-b[0]), fmax(fabs(a[1]-b[1]), fabs(a[2]-b[2])));
}

// Line-integral field of loop lp (polyline P, fan apex c) at p, no periodic image
static void lineField(const loop& lp, const vector<V>& P, const V& c, const V& p, double u[3])
{
    vector<double> X, Y, Z;
    for (const V& q : P) { X.push_back(q[0]); Y.push_back(q[1]); Z.push_back(q[2]); }
    const V m = lp.getM();
    const double bm = lp.getBurgers()[0];
    u[0] = u[1] = u[2] = 0.;
    barnettLoop<double>(p[0], p[1], p[2], X.data(), Y.data(), Z.data(), (int)P.size(),
                        c[0], c[1], c[2], bm*m[0], bm*m[1], bm*m[2], u[0], u[1], u[2]);
}

// Kroupa kernel integrated over the fan triangles (c, P_k, P_k+1), each cut into
// sub^2 triangles with the 7-point degree-5 rule (Dunavant)
static void kroupaQuadrature(const loop& lp, const vector<V>& P, const V& c, const V& p,
                             int sub, double u[3])
{
    const V m = lp.getM(), n = lp.getN();
    const double bm = lp.getBurgers()[0];
    LoopPhys L = { bm, bm, bm, m[0], m[1], m[2], n[0], n[1], n[2], 0. };
    const double k1 = kroupaK1(), bn = dot(m, n);
    const double W[7] = { 0.225, 0.1323941527, 0.1323941527, 0.1323941527,
                          0.1259391805, 0.1259391805, 0.1259391805 };
    const double a1 = 0.0597158717, b1 = 0.4701420641, a2 = 0.7974269853, b2 = 0.1012865073;
    const double B[7][3] = { {1/3., 1/3., 1/3.}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1},
                             {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2} };
    u[0] = u[1] = u[2] = 0.;
    const int N = (int)P.size();
    for (int t = 0; t < N; t++) {
        const V& A = P[t];
        const V& Bv = P[(t + 1) % N];
        const V e1 = { (A[0]-c[0])/sub, (A[1]-c[1])/sub, (A[2]-c[2])/sub };
        const V e2 = { (Bv[0]-c[0])/sub, (Bv[1]-c[1])/sub, (Bv[2]-c[2])/sub };
        const V cc = cross(e1, e2);
        const double area = 0.5 * sqrt(dot(cc, cc));
        auto pt = [&](double i, double j) {
            return V{ c[0] + i*e1[0] + j*e2[0], c[1] + i*e1[1] + j*e2[1], c[2] + i*e1[2] + j*e2[2] };
        };
        for (int i = 0; i < sub; i++)
            for (int j = 0; j < sub - i; j++)
                for (int up = 0; up < 2; up++) {
                    if (up && i + j >= sub - 1) continue;
                    const V q0 = up ? pt(i+1, j) : pt(i, j);
                    const V q1 = up ? pt(i+1, j+1) : pt(i+1, j);
                    const V q2 = pt(i, j+1);
                    for (int g = 0; g < 7; g++) {
                        double x[3];
                        for (int d = 0; d < 3; d++) x[d] = B[g][0]*q0[d] + B[g][1]*q1[d] + B[g][2]*q2[d];
                        L.dS = W[g] * area;
                        kroupaPair(p[0]-x[0], p[1]-x[1], p[2]-x[2], L, k1, bn, u[0], u[1], u[2]);
                    }
                }
    }
}

// Test points around loop lp: heights |s| in [1, 3] Ang (both sides), in-plane
// offsets up to 1.3 x size, deterministic
static vector<V> testPoints(const loop& lp, double size, int count)
{
    const V c = lp.getCenter(), n = lp.getN(), mm = lp.getM();
    V e1 = { mm[0] - dot(mm, n)*n[0], mm[1] - dot(mm, n)*n[1], mm[2] - dot(mm, n)*n[2] };
    if (dot(e1, e1) < 1e-12) e1 = cross(n, V{ 1., 0., 0. });
    if (dot(e1, e1) < 1e-12) e1 = cross(n, V{ 0., 1., 0. });
    const double l = sqrt(dot(e1, e1));
    for (double& x : e1) x /= l;
    const V e2 = cross(n, e1);
    vector<V> pts;
    unsigned seed = 2024;
    auto rnd = [&]() { seed = seed * 1103515245u + 12345u; return (seed >> 8) / double(1u << 24); };
    for (int k = 0; k < count; k++) {
        const double a = (2.*rnd() - 1.) * 1.3 * size, b = (2.*rnd() - 1.) * 1.3 * size;
        const double h = (1. + 2.*rnd()) * (rnd() < 0.5 ? -1. : 1.);
        pts.push_back({ c[0] + a*e1[0] + b*e2[0] + h*n[0], c[1] + a*e1[1] + b*e2[1] + h*n[1],
                        c[2] + a*e1[2] + b*e2[2] + h*n[2] });
    }
    return pts;
}

// 1 + 2 + 3 + 6 for one loop
static void checkLoop(const char* label, loop lp, double size, int sub)
{
    const double eps = 0.01;
    const vector<V> P = lp.polyline(eps);
    const V c = lp.getCenter(), n = lp.getN(), m = lp.getM();
    const double bm = lp.getBurgers()[0];

    // 6. polyline: CCW about n, area close to the exact contour area
    double area = 0.;
    for (size_t k = 0; k < P.size(); k++) {
        const V& p = P[k];
        const V& q = P[(k + 1) % P.size()];
        area += 0.5 * dot(cross(V{ p[0]-c[0], p[1]-c[1], p[2]-c[2] }, V{ q[0]-c[0], q[1]-c[1], q[2]-c[2] }), n);
    }
    assert(area > 0.);
    const double exact = lp.getContourArea();
    assert(fabs(area - exact) < 2e-3 * exact);

    // 1. against the Kroupa quadrature
    double err = 0., umax = 0.;
    for (const V& p : testPoints(lp, size, 8)) {
        double ul[3], uq[3];
        lineField(lp, P, c, p, ul);
        kroupaQuadrature(lp, P, c, p, sub, uq);
        const double z[3] = { 0., 0., 0. };
        err  = fmax(err, maxAbs3(ul, uq));
        umax = fmax(umax, maxAbs3(uq, z));
    }
    assert(err < 1e-6);

    // 2. jump across the cut: -b inside, 0 outside (well clear of the contour)
    const double h = 1e-7;
    const V pin  = { c[0] + 0.3*size*m[0]*(1. - fabs(dot(m, n))), c[1], c[2] };
    double up[3], um[3], jumpErr = 0.;
    for (int side = 0; side < 2; side++) {
        // inside: the centre; outside: 3 x size away in the plane
        V q = c;
        if (side == 1) {
            V e = cross(n, V{ 0.3, 0.5, 0.8 });
            const double l = sqrt(dot(e, e));
            for (int i = 0; i < 3; i++) q[i] += 3. * size * e[i] / l;
        }
        lineField(lp, P, c, V{ q[0] + h*n[0], q[1] + h*n[1], q[2] + h*n[2] }, up);
        lineField(lp, P, c, V{ q[0] - h*n[0], q[1] - h*n[1], q[2] - h*n[2] }, um);
        for (int i = 0; i < 3; i++) {
            const double expect = (side == 0) ? -bm * m[i] : 0.;
            jumpErr = fmax(jumpErr, fabs(up[i] - um[i] - expect));
        }
    }
    (void)pin;
    assert(jumpErr < 1e-5);

    // 3. invariance: apex moved inside the contour, starting vertex rotated
    double inv = 0.;
    vector<V> Pr(P.begin() + P.size() / 3, P.end());
    Pr.insert(Pr.end(), P.begin(), P.begin() + P.size() / 3);
    V c2 = c;
    {
        V e = cross(n, V{ 0.7, -0.2, 0.4 });
        const double l = sqrt(dot(e, e));
        for (int i = 0; i < 3; i++) c2[i] += 0.25 * size * e[i] / l;
    }
    for (const V& p : testPoints(lp, size, 6)) {
        double u1[3], u2[3];
        lineField(lp, P, c, p, u1);
        lineField(lp, Pr, c2, p, u2);
        inv = fmax(inv, maxAbs3(u1, u2));
    }
    assert(inv < 1e-11);

    cout << "  " << label << ": " << P.size() << " segments, area " << area << " (exact " << exact
         << "), vs quadrature " << err << " Ang (|u| up to " << umax << "), jump error " << jumpErr
         << ", invariance " << inv << "\n";
}

// 4. Kroupa loopon sum -> line integral as da -> 0 (sharp hexagon: same contour)
static void checkLooponConvergence()
{
    const V m = { -1, -1, -1 }, n = { 1, 1, 1 }, c = { 0.2, -0.1, 0.3 };
    double prev = 1e9;
    cout << "  loopon sum vs line integral, prismatic hexagon R = 10 Ang (points >= 1.5 Ang off the plane):\n";
    for (double da : { 0.8, 0.4, 0.2, 0.1 }) {
        loop lp(da, regularPolygon(6, 10., {}, m, n), 0., c, m, n, A0, ONE);
        lp.generateSurface();
        lp.generateLoopons();
        const vector<V> P = lp.polyline(0.01);
        const V mh = lp.getM(), nh = lp.getN();
        const double bm = lp.getBurgers()[0];
        const LoopPhys L = { bm, bm, bm, mh[0], mh[1], mh[2], nh[0], nh[1], nh[2], lp.getSurface() };
        double err = 0.;
        for (const V& p0 : testPoints(lp, 10., 10)) {
            // keep points at least 1.5 Ang from the plane
            const double s = dot(V{ p0[0]-c[0], p0[1]-c[1], p0[2]-c[2] }, nh);
            const double shift = (fabs(s) < 1.5) ? (s > 0 ? 1.5 - s : -1.5 - s) : 0.;
            const V p = { p0[0] + shift*nh[0], p0[1] + shift*nh[1], p0[2] + shift*nh[2] };
            double ul[3], uk[3] = { 0., 0., 0. };
            lineField(lp, P, lp.getCenter(), p, ul);
            for (int k = 0; k < lp.getNumberOfLoopons(); k++) {
                const V q = lp.getOneLoopon(k);
                kroupaPair(p[0]-q[0], p[1]-q[1], p[2]-q[2], L, kroupaK1(), dot(mh, nh), uk[0], uk[1], uk[2]);
            }
            err = fmax(err, maxAbs3(ul, uk));
        }
        cout << "    da = " << da << " Ang: max |u_loopons - u_line| = " << err << " Ang\n";
        assert(err < prev);
        prev = err;
    }
    assert(prev < 2e-3);
}

// 5. Far field: one loopon of the loop's area at its centre
static void checkFarField()
{
    const V m = { 1, -1, 1 }, n = { 1, 1, 0 }, c = { 0., 0., 0. };
    loop lp(0.1, 2.0, c, m, n, A0, ONE);
    const vector<V> P = lp.polyline(0.001);
    const V mh = lp.getM(), nh = lp.getN();
    const double bm = lp.getBurgers()[0];
    const LoopPhys L = { bm, bm, bm, mh[0], mh[1], mh[2], nh[0], nh[1], nh[2], lp.getContourArea() };
    const V p = { 60., -35., 80. };
    double ul[3], uk[3] = { 0., 0., 0. };
    lineField(lp, P, c, p, ul);
    kroupaPair(p[0], p[1], p[2], L, kroupaK1(), dot(mh, nh), uk[0], uk[1], uk[2]);
    const double z[3] = { 0., 0., 0. };
    const double rel = maxAbs3(ul, uk) / maxAbs3(uk, z);
    cout << "  far field (R = 2 Ang loop at 105 Ang): relative difference to one loopon " << rel << "\n";
    assert(rel < 1e-3);
}

int main()
{
    cout << "test_barnett:\n";
    const V c = { 0.3, -0.2, 0.1 };
    // glide ({110}), prismatic ({111}), mixed (1/2<111> on {110}); several shapes
    checkLoop("circle, glide        ", loop(0.4, 8., c, { 1, -1, 1 }, { 1, 1, 0 }, A0, ONE), 8., 30);
    checkLoop("circle, prismatic    ", loop(0.4, 8., c, { -1, -1, -1 }, { 1, 1, 1 }, A0, ONE), 8., 30);
    checkLoop("rounded rect, mixed  ", loop(0.4, 9., 6., 2.5, c, { 1, 1, 1 }, { 1, -1, 0 }, A0, ONE), 9., 40);
    checkLoop("parallelogram, glide ", loop(0.4, 7., 5., { 1, 1, 1 }, { 1, 1, -1 }, c, { 1, 1, 1 },
                                            { 1, -1, 0 }, A0, ONE), 9., 60);
    checkLoop("hexagon, interstitial", loop(0.4, regularPolygon(6, 8., { 1, -1, 0 }, { -1, -1, -1 }, { 1, 1, 1 }),
                                            0., c, { -1, -1, -1 }, { 1, 1, 1 }, A0, ONE, { 1, -1, 0 }), 8., 80);
    checkLoop("hexagon rounded, vac.", loop(0.4, regularPolygon(6, 8., {}, { 1, 1, 1 }, { 1, 1, 1 }),
                                            2., c, { 1, 1, 1 }, { 1, 1, 1 }, A0, ONE), 8., 40);
    checkLoop("a[001] square, prism.", loop(0.4, 6., 6., 0., c, { 0, 0, 2 }, { 0, 0, 1 }, A0, ONE), 6., 60);
    checkLooponConvergence();
    checkFarField();
    cout << "test_barnett: OK\n";
    return 0;
}
