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

int main()
{
    cout << "test_polycrystal:\n";
    oneGrainEqualsOriented();
    randomGrains(Lattice::BCC);
    randomGrains(Lattice::FCC);

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
