// =============================================================================
// Unit test for the periodic Voronoi polycrystal (polycrystal.hpp):
//   1. One grain with a commensurate orientation, seeded on the box corner,
//      gives exactly the sites of generateOrientedLattice (half-open cells:
//      no site kept twice, none missing, no GB removal).
//   2. Random grains (BCC, FCC): every atom lies in the Voronoi cell of its
//      grain (grainOf), no pair is closer than r_min, the atom count is within
//      a few % below the perfect crystal, and the build is deterministic.
//   3. Rotations: random and Euler orientations are proper rotations; Euler
//      angles round-trip.
//   4. PortableRandom reproduces the standard mt19937_64 sequence.
//   5. Loop clipping (cellHalfSpaces, clipPolygon): a polygon inside a cell is
//      unchanged; a polygon larger than the box becomes the cell's cross-section
//      (vertices on or inside the cell, area = brute-force grid count, convex,
//      same orientation); a polygon away from the cell is empty.
//
//   make test
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../polycrystal.hpp"

using namespace std;

static void checkRotation(const Orientation& o)
{
    double det = 0.;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double dot = 0.;
            for (int k = 0; k < 3; k++) dot += o.R[i][k] * o.R[j][k];
            assert(fabs(dot - (i == j ? 1. : 0.)) < 1e-12);
        }
    det = o.R[0][0]*(o.R[1][1]*o.R[2][2] - o.R[1][2]*o.R[2][1])
        - o.R[0][1]*(o.R[1][0]*o.R[2][2] - o.R[1][2]*o.R[2][0])
        + o.R[0][2]*(o.R[1][0]*o.R[2][1] - o.R[1][1]*o.R[2][0]);
    assert(fabs(det - 1.) < 1e-12);
}

static void oneGrainEqualsOriented()
{
    const double a0 = 3.6;
    Orientation o;
    const long d[3][3] = { {1, -1, 0}, {1, 1, -2}, {1, 1, 1} };
    o.active = true;
    for (int a = 0; a < 3; a++) {
        const double n = sqrt(double(d[a][0]*d[a][0] + d[a][1]*d[a][1] + d[a][2]*d[a][2]));
        for (int k = 0; k < 3; k++) { o.dir[a][k] = d[a][k]; o.R[a][k] = d[a][k] / n; }
    }
    PolycrystalSpec spec;
    spec.active = true;
    const int np[3] = { 6, 4, 3 };
    for (int a = 0; a < 3; a++) spec.box[a] = np[a] * latticePeriod(Lattice::FCC, o.dir[a], a0);
    Grain g;
    for (int k = 0; k < 3; k++) g.seed[k] = -0.5 * spec.box[k];   // lattice site on the box corner
    g.o = o;
    spec.list = { g };
    spec.hasCenter = { 1 };
    spec.hasOrientation = { 1 };
    const PolycrystalResult r = buildPolycrystal(Lattice::FCC, a0, spec, makeGrains(spec));
    const vector<vector<double>> ref = generateOrientedLattice(Lattice::FCC, o, spec.box, a0);
    assert(r.removedGB == 0);
    assert(r.x.size() == ref.size());
    auto key = [](double x, double y, double z) { return vector<long>{ lround(x * 1e5), lround(y * 1e5), lround(z * 1e5) }; };
    vector<vector<long>> A, B;
    for (size_t i = 0; i < ref.size(); i++) {
        A.push_back(key(r.x[i], r.y[i], r.z[i]));
        B.push_back(key(ref[i][0], ref[i][1], ref[i][2]));
    }
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    assert(A == B);
    cout << "  1 grain, FCC X=[1-10] Y=[11-2] Z=[111]: " << ref.size() << " atoms = oriented lattice\n";
}

static void randomGrains(Lattice lat)
{
    const double a0 = 3.4;
    PolycrystalSpec spec;
    spec.active = true;
    spec.box[0] = 44.; spec.box[1] = 38.; spec.box[2] = 41.;
    spec.nGrains = 6;
    spec.seed = 11;
    const vector<Grain> grains = makeGrains(spec);
    const PolycrystalResult r = buildPolycrystal(lat, a0, spec, grains);
    const long N = (long)r.x.size();
    const double* L = spec.box;
    // 1. own Voronoi cell (ties are measure-zero for random seeds)
    for (long a = 0; a < N; a++) {
        const double x[3] = { r.x[a], r.y[a], r.z[a] };
        assert(grainOf(grains, L, x) == r.grain[a]);
        for (int k = 0; k < 3; k++) assert(x[k] >= -0.5 * L[k] && x[k] < 0.5 * L[k]);
    }
    // 2. no pair closer than r_min (brute force, periodic)
    const double rmin = spec.gbMinDistance * nearestNeighbourDistance(lat, a0);
    double dmin = 1e9;
    for (long a = 0; a < N; a++)
        for (long b = a + 1; b < N; b++) {
            double s2 = 0.;
            const double e[3] = { r.x[b] - r.x[a], r.y[b] - r.y[a], r.z[b] - r.z[a] };
            for (int k = 0; k < 3; k++) { const double v = e[k] - L[k] * floor(e[k] / L[k] + 0.5); s2 += v * v; }
            dmin = min(dmin, s2);
        }
    dmin = sqrt(dmin);
    assert(dmin >= rmin - 1e-9);
    // 3. before GB removal the Voronoi fill holds the perfect-crystal density (each
    //    site of each grain's lattice counted once); the removal itself is large
    //    here (grains of ~20 Ang are mostly boundary) but bounded
    const double perfect = atomsPerCell(lat) * L[0] * L[1] * L[2] / (a0 * a0 * a0);
    assert(fabs((N + r.removedGB) - perfect) < 0.01 * perfect);
    assert(r.removedGB < 0.2 * perfect);
    long sum = 0;
    for (const auto& g : r.grains) sum += g.atoms;
    assert(sum == N);
    // 4. deterministic
    const PolycrystalResult r2 = buildPolycrystal(lat, a0, spec, makeGrains(spec));
    assert(r2.x == r.x && r2.grain == r.grain);
    for (const auto& g : grains) checkRotation(g.o);
    cout << "  " << latticeName(lat) << " 6 random grains: " << N << " atoms (perfect " << lround(perfect)
         << ", " << r.removedGB << " removed at GBs), min distance " << dmin << " >= r_min " << rmin << "\n";
}

static vector<vector<double>> square(const double c[3], const double u[3], const double v[3], double h)
{
    vector<vector<double>> P;
    const double su[4] = { 1, -1, -1, 1 }, sv[4] = { 1, 1, -1, -1 };
    for (int q = 0; q < 4; q++)
        P.push_back({ c[0] + h * (su[q]*u[0] + sv[q]*v[0]), c[1] + h * (su[q]*u[1] + sv[q]*v[1]),
                      c[2] + h * (su[q]*u[2] + sv[q]*v[2]) });
    return P;
}

static void clipping()
{
    PolycrystalSpec spec;
    spec.active = true;
    spec.box[0] = 60.; spec.box[1] = 50.; spec.box[2] = 55.;
    spec.nGrains = 7;
    spec.seed = 5;
    const vector<Grain> grains = makeGrains(spec);
    const double u[3] = { 0.6, 0.8, 0. }, v[3] = { -0.48, 0.36, 0.8 };   // orthonormal
    const double nrm[3] = { u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] };
    for (int k = 0; k < (int)grains.size(); k++) {
        const vector<HalfSpace> hs = cellHalfSpaces(grains, spec.box, k);
        const double* s = grains[k].seed;
        // inside: a 1 Ang square around the seed is unchanged
        const auto small = square(s, u, v, 0.5);
        const auto ks = clipPolygon(small, hs);
        assert(ks.size() == 4 && fabs(polygonArea(ks) - 1.) < 1e-12);
        // larger than the box: the cell cross-section
        const double H = 100.;
        const auto big = clipPolygon(square(s, u, v, H), hs);
        assert(big.size() >= 3);
        for (const auto& x : big) { const double xx[3] = { x[0], x[1], x[2] }; assert(insideCell(hs, xx, 1e-6)); }
        double ori = 0.;                       // orientation kept (normal along +nrm), convex
        for (size_t q = 0; q < big.size(); q++) {
            const auto& A = big[q]; const auto& B = big[(q + 1) % big.size()]; const auto& C = big[(q + 2) % big.size()];
            const double e1[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] }, e2[3] = { C[0]-B[0], C[1]-B[1], C[2]-B[2] };
            const double t = (e1[1]*e2[2] - e1[2]*e2[1]) * nrm[0] + (e1[2]*e2[0] - e1[0]*e2[2]) * nrm[1]
                           + (e1[0]*e2[1] - e1[1]*e2[0]) * nrm[2];
            assert(t > 0.);
            ori += t;
        }
        assert(ori > 0.);
        const int M = 800;                    // grid count of the cross-section
        const double step = 2. * H / M;
        long in = 0;
        for (int a = 0; a < M; a++)
            for (int b = 0; b < M; b++) {
                const double ca = -H + (a + 0.5) * step, cb = -H + (b + 0.5) * step;
                const double x[3] = { s[0] + ca*u[0] + cb*v[0], s[1] + ca*u[1] + cb*v[1], s[2] + ca*u[2] + cb*v[2] };
                in += insideCell(hs, x, 0.);
            }
        const double gridArea = in * step * step, area = polygonArea(big);
        assert(fabs(gridArea - area) < 0.01 * area + 4. * step * sqrt(area));
        // away from the cell: one box length out along the plane, beyond the cell
        const double far[3] = { s[0] + 3. * spec.box[0] * u[0], s[1] + 3. * spec.box[0] * u[1], s[2] + 3. * spec.box[0] * u[2] };
        assert(clipPolygon(square(far, u, v, 5.), hs).empty());
        if (k == 0)
            cout << "  clipping, grain 1: cross-section " << big.size() << " vertices, area " << area
                 << " Ang^2 (grid " << gridArea << ")\n";
    }
    cout << "  clipping: 7 grains, inside / cross-section / outside OK\n";
}

int main()
{
    cout << "test_polycrystal:\n";
    oneGrainEqualsOriented();
    randomGrains(Lattice::BCC);
    randomGrains(Lattice::FCC);
    clipping();

    // Euler round trip
    for (const double* e : { (const double[]){ 30., 40., 50. }, (const double[]){ 200., 120., 10. },
                             (const double[]){ 359., 90., 270. } }) {
        const Orientation o = orientationFromEuler(e[0], e[1], e[2]);
        checkRotation(o);
        double p1, P, p2;
        eulerFromOrientation(o, p1, P, p2);
        assert(fabs(p1 - e[0]) < 1e-9 && fabs(P - e[1]) < 1e-9 && fabs(p2 - e[2]) < 1e-9);
    }
    // Euler (0, 0, 0) is the cube orientation
    const Orientation id = orientationFromEuler(0., 0., 0.);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) assert(fabs(id.R[i][j] - (i == j)) < 1e-15);

    // mt19937_64: the 10000th output for the default seed 5489 is fixed by the standard
    std::mt19937_64 eng;
    eng.discard(9999);
    assert(eng() == 9981545732273789042ULL);
    PortableRandom pr(5489);
    assert(pr.uniform() == (double)(14514284786278117030ULL >> 11) * (1.0 / 9007199254740992.0));
    cout << "test_polycrystal: OK\n";
    return 0;
}
