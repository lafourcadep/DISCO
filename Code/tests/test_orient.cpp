// =============================================================================
// Unit test for oriented crystals (lattice.hpp, snapToMidPlaneOriented):
//   1. generateOrientedLattice: atom count = atoms per cell x V / a0^3, and every
//      atom has exactly the perfect first-neighbour shell (12 FCC, 8 BCC) at the
//      nearest-neighbour distance, periodic boundaries included: no duplicated or
//      missing site on the box faces, so the box is commensurate and periodic.
//   2. snapToMidPlaneOriented: every atom of the generated lattice lies at a
//      half-integer number of plane spacings from the snapped centre, for
//      crystal-frame normals of several families, and the shift is along n.
//   3. rotateToLab / rotateToCrystal are inverse to each other.
//
//   make test
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../loop.hpp"

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

static void check(Lattice lat, const long d[3][3])
{
    const double a0 = 3.6;
    const Orientation o = makeOrientation(d);
    double L[3];
    for (int a = 0; a < 3; a++) {
        const double p = latticePeriod(lat, o.dir[a], a0);
        L[a] = ceil(11. / p) * p;                     // >= 11 Ang, whole periods
    }
    const vector<vector<double>> at = generateOrientedLattice(lat, o, L, a0);
    const double expected = atomsPerCell(lat) * L[0] * L[1] * L[2] / (a0 * a0 * a0);
    assert(fabs(expected - lround(expected)) < 1e-6);
    assert((long)at.size() == lround(expected));

    // 1. first-neighbour shell, periodic
    const double nn = (lat == Lattice::FCC) ? a0 / sqrt(2.) : a0 * sqrt(3.) / 2.;
    const int zWant = (lat == Lattice::FCC) ? 12 : 8;
    const size_t N = at.size();
    for (size_t i = 0; i < N; i++) {
        int z = 0;
        for (size_t j = 0; j < N; j++) {
            if (i == j) continue;
            double r2 = 0.;
            for (int k = 0; k < 3; k++) {
                double r = at[j][k] - at[i][k];
                r -= L[k] * floor(r / L[k] + 0.5);
                r2 += r * r;
            }
            const double r = sqrt(r2);
            assert(r > nn * (1. - 1e-6));             // no duplicate or too-close site
            if (r < nn * (1. + 1e-6)) z++;
        }
        assert(z == zWant);
    }

    // 2. oriented mid-plane snap, brute force over the generated sites
    const double origin[3] = { -L[0] / 2, -L[1] / 2, -L[2] / 2 };
    for (const vector<double>& n : vector<vector<double>>{ {1, 1, 1}, {1, -1, 0}, {1, 1, -2}, {0, 0, 1}, {1, 2, -3} }) {
        const vector<double> c = { 0.37, -1.21, 0.83 };
        double sp;
        const vector<double> s = snapToMidPlaneOriented(c, n, o, origin, a0, sp, lat);
        const double ln = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        const vector<double> nh = rotateToLab(o, { n[0]/ln, n[1]/ln, n[2]/ln });
        const double sMid = s[0]*nh[0] + s[1]*nh[1] + s[2]*nh[2];
        const double sC   = c[0]*nh[0] + c[1]*nh[1] + c[2]*nh[2];
        for (int k = 0; k < 3; k++) assert(fabs((s[k] - c[k]) - (sMid - sC) * nh[k]) < 1e-9);
        assert(fabs(sMid - sC) <= sp + 1e-9);
        double nearest = 1e9;
        for (const auto& p : at) {
            const double t = ((p[0]*nh[0] + p[1]*nh[1] + p[2]*nh[2]) - sMid) / sp - 0.5;
            assert(fabs(t - lround(t)) < 1e-6);       // every atom on a plane, centre between two
            nearest = min(nearest, fabs((p[0]*nh[0] + p[1]*nh[1] + p[2]*nh[2]) - sMid));
        }
        assert(fabs(nearest - sp / 2) < 1e-6);
    }

    // 3. rotations
    const vector<double> v = { 0.3, -1.7, 2.2 };
    const vector<double> w = rotateToCrystal(o, rotateToLab(o, v));
    for (int k = 0; k < 3; k++) assert(fabs(w[k] - v[k]) < 1e-12);

    cout << "  " << latticeName(lat) << " X=[" << d[0][0] << " " << d[0][1] << " " << d[0][2] << "] Y=["
         << d[1][0] << " " << d[1][1] << " " << d[1][2] << "] Z=[" << d[2][0] << " " << d[2][1] << " " << d[2][2]
         << "]: box " << L[0] << " x " << L[1] << " x " << L[2] << ", " << N << " atoms, shells and planes OK\n";
}

int main()
{
    cout << "test_orient:\n";
    const long orients[][3][3] = {
        { {1, 0, 0},  {0, 1, 0},  {0, 0, 1} },
        { {1, -1, 0}, {1, 1, -2}, {1, 1, 1} },        // FCC dislocation MD standard
        { {1, 1, 1},  {-1, 1, 0}, {-1, -1, 2} },      // BCC screw standard
        { {1, 1, 0},  {-1, 1, 0}, {0, 0, 1} },
        { {1, 2, 0},  {-2, 1, 0}, {0, 0, 1} },
    };
    for (Lattice lat : { Lattice::BCC, Lattice::FCC })
        for (const auto& d : orients) check(lat, d);
    // Periods along common directions
    const double a0 = 3.6;
    const long d110[3] = { 1, 1, 0 }, d111[3] = { 1, 1, 1 }, d112[3] = { 1, 1, -2 }, d100[3] = { 2, 0, 0 };
    assert(fabs(latticePeriod(Lattice::FCC, d110, a0) - a0 / sqrt(2.)) < 1e-12);
    assert(fabs(latticePeriod(Lattice::FCC, d111, a0) - a0 * sqrt(3.)) < 1e-12);
    assert(fabs(latticePeriod(Lattice::FCC, d112, a0) - a0 * sqrt(6.) / 2) < 1e-12);
    assert(fabs(latticePeriod(Lattice::BCC, d111, a0) - a0 * sqrt(3.) / 2) < 1e-12);
    assert(fabs(latticePeriod(Lattice::BCC, d110, a0) - a0 * sqrt(2.)) < 1e-12);
    assert(fabs(latticePeriod(Lattice::BCC, d100, a0) - a0) < 1e-12);
    cout << "test_orient: OK\n";
    return 0;
}
