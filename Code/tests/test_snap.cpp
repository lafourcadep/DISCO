// =============================================================================
// Unit test for snapToMidPlane: the snapped centre must sit exactly halfway
// between two atomic planes of the (hkl) family, checked by brute force over
// the BCC sites of generateBCC(), with in-plane components unchanged.
//
//   make test
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include "../loop.hpp"

using namespace std;

static void check(vector<double> n, int nx, int ny, int nz, vector<double> c)
{
    const double a0 = 3.3024;
    double d;
    vector<double> s = snapToMidPlane(c, n, a0, nx, ny, nz, d);

    double nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    double nh[3] = { n[0]/nn, n[1]/nn, n[2]/nn };
    double ps = s[0]*nh[0] + s[1]*nh[1] + s[2]*nh[2];

    // Shift is purely along n
    double pc = c[0]*nh[0] + c[1]*nh[1] + c[2]*nh[2];
    for (int i = 0; i < 3; i++)
        assert(fabs((s[i] - c[i]) - (ps - pc) * nh[i]) < 1e-9);

    // Nearest atomic planes below and above the snapped centre (same sites as generateBCC)
    double below = -1e30, above = 1e30;
    const double hx = nx*a0/2., hy = ny*a0/2., hz = nz*a0/2.;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                for (double f : {0.0, 0.5}) {
                    double p = (i*a0 + f*a0 - hx)*nh[0] + (j*a0 + f*a0 - hy)*nh[1]
                             + (k*a0 + f*a0 - hz)*nh[2];
                    if (p <= ps) below = max(below, p);
                    else         above = min(above, p);
                }
    assert(fabs((above - below) - d) < 1e-9);      // reported spacing is the real one
    assert(fabs((ps - below) - d/2) < 1e-9);       // exactly on the mid-plane
    cout << "  n = (" << n[0] << " " << n[1] << " " << n[2] << "), N = "
         << nx << "x" << ny << "x" << nz << ": d = " << d << " Ang, mid-plane OK\n";
}

int main()
{
    const vector<double> c = { 1.3, -2.7, 0.9 };       // arbitrary centre (Ang)
    for (auto dims : vector<vector<int>>{ {12, 12, 12}, {11, 13, 9} })
        for (auto n : vector<vector<double>>{ {1, 1, 0}, {1, -1, 2}, {1, 2, -3}, {1, 0, 0}, {2, 2, 0} })
            check(n, dims[0], dims[1], dims[2], c);

    // A centre exactly on an atomic plane goes to the mid-plane on the +n side
    double d;
    vector<double> onPlane = { 0.0, 0.0, 0.0 };         // BCC site for even N
    vector<double> s = snapToMidPlane(onPlane, {1, 1, 0}, 3.3024, 12, 12, 12, d);
    assert(fabs((s[0] + s[1]) / sqrt(2.) - d / 2) < 1e-9);

    cout << "All snap tests passed.\n";
    return 0;
}
