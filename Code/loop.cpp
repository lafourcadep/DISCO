#ifndef LOOP_CPP
#define LOOP_CPP

#include "loop.hpp"

using namespace std;

// =============================================================================
// Shared initialisation helper (avoids code duplication between constructors)
// =============================================================================
static void initCommon(
    vector<double>& _m, const vector<double>& m,
    vector<double>& _n, const vector<double>& n,
    vector<double>& _p,
    vector<double>& _burgers,
    const vector<double>& _compress,
    double _a0)
{
    double normM = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
    double normN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    for (int i = 0; i < 3; i++) { _m[i] = m[i] / normM; }
    for (int i = 0; i < 3; i++) { _n[i] = n[i] / normN; }

    // Burgers vector magnitude per Cartesian direction: |b| = a0 |m| / 2
    for (int j = 0; j < 3; j++)
        _burgers[j] = _compress[j] * _a0 * 0.5 * normM;

    // p = m × n  (third basis vector, in-plane, perpendicular to m)
    _p[0] = _m[1]*_n[2] - _m[2]*_n[1];
    _p[1] = _m[2]*_n[0] - _m[0]*_n[2];
    _p[2] = _m[0]*_n[1] - _m[1]*_n[0];
}

// =============================================================================
// Circular constructor
// =============================================================================
loop::loop(double da, double radius,
           vector<double> center, vector<double> m, vector<double> n,
           double latticeParameter, vector<double> compress)
  : _shape(LoopShape::CIRCULAR),
    _radius(radius),
    _half_m(0.), _half_p(0.), _r_corner(0.),
    _center(center), _m(m), _n(n),
    _burgers(3, 0.), _p(3, 0.),
    _compress(compress), _a0(latticeParameter), _da(da)
{
    initCommon(_m, m, _n, n, _p, _burgers, _compress, _a0);
    _u = _m;
    _v = _p;
}

// =============================================================================
// Rectangular constructor
// Without `side`: sides parallel to m = pure screw, parallel to p = pure edge.
// With `side` (in-plane, checked by input.cpp): the half_m axis is along side
// and the half_p axis along side × n; b stays along m (mixed-character sides).
// =============================================================================
loop::loop(double da, double half_m, double half_p, double corner_r,
           vector<double> center, vector<double> m, vector<double> n,
           double latticeParameter, vector<double> compress,
           vector<double> side)
  : _shape(LoopShape::RECTANGULAR),
    _radius(0.),
    _half_m(half_m), _half_p(half_p),
    _center(center), _m(m), _n(n),
    _burgers(3, 0.), _p(3, 0.),
    _compress(compress), _a0(latticeParameter), _da(da)
{
    initCommon(_m, m, _n, n, _p, _burgers, _compress, _a0);
    if (side.empty()) {
        _u = _m;
        _v = _p;
    } else {
        double norm = sqrt(side[0]*side[0] + side[1]*side[1] + side[2]*side[2]);
        _u = { side[0] / norm, side[1] / norm, side[2] / norm };
        _v = { _u[1]*_n[2] - _u[2]*_n[1],
               _u[2]*_n[0] - _u[0]*_n[2],
               _u[0]*_n[1] - _u[1]*_n[0] };
    }

    // Clamp corner radius to the smaller half-dimension
    double maxCorner = min(_half_m, _half_p);
    if (corner_r > maxCorner) {
        cerr << "WARNING: corner_radius " << corner_r
             << " clamped to min(Half_m, Half_p) = " << maxCorner << "\n";
        corner_r = maxCorner;
    }
    _r_corner = corner_r;
}

// =============================================================================
// Parallelogram constructor
// =============================================================================
loop::loop(double da, double half1, double half2,
           vector<double> dir1, vector<double> dir2,
           vector<double> center, vector<double> m, vector<double> n,
           double latticeParameter, vector<double> compress)
  : _shape(LoopShape::PARALLELOGRAM),
    _radius(0.),
    _half_m(half1), _half_p(half2), _r_corner(0.),
    _center(center), _m(m), _n(n),
    _burgers(3, 0.), _p(3, 0.),
    _compress(compress), _a0(latticeParameter), _da(da)
{
    initCommon(_m, m, _n, n, _p, _burgers, _compress, _a0);
    const double l1 = sqrt(dir1[0]*dir1[0] + dir1[1]*dir1[1] + dir1[2]*dir1[2]);
    const double l2 = sqrt(dir2[0]*dir2[0] + dir2[1]*dir2[1] + dir2[2]*dir2[2]);
    _u = { dir1[0] / l1, dir1[1] / l1, dir1[2] / l1 };
    _v = { dir2[0] / l2, dir2[1] / l2, dir2[2] / l2 };
}

// =============================================================================
// generateLoopons
// All shapes use the same lattice u = i*da, v = j*da (in-plane, Angstrom)
// along the axes _u, _v (oblique for a parallelogram), centred on the loop
// centre; a node is kept if it lies inside the loop. A parallelogram uses the
// rectangle test in its oblique coordinates, with no corner rounding.
// =============================================================================
void loop::generateLoopons()
{
    vector<double> singleLoopon(3);

    // Grid points exactly on the contour: excluded for circles, included for
    // rectangles (as in the original code). The tolerance makes this robust to
    // rounding once sizes are given in Angstrom.
    const double eps = 1e-9;
    const double inner_m = max(_half_m - _r_corner, 0.0);
    const double inner_p = max(_half_p - _r_corner, 0.0);
    const double rc2     = _r_corner * _r_corner;

    const double ext_m = (_shape == LoopShape::CIRCULAR) ? _radius : _half_m;
    const double ext_p = (_shape == LoopShape::CIRCULAR) ? _radius : _half_p;
    const int steps_m = (int)ceil(ext_m / _da) + 1;
    const int steps_p = (int)ceil(ext_p / _da) + 1;

    for (int xi = -steps_m; xi <= steps_m; xi++) {
        for (int yi = -steps_p; yi <= steps_p; yi++) {
            double u = xi * _da;
            double v = yi * _da;
            bool inside;
            if (_shape == LoopShape::CIRCULAR) {
                inside = sqrt(u*u + v*v) < _radius - eps;
            } else {
                // Rounded rectangle:
                //   max(|u| - (half_m - r_corner), 0)² + max(|v| - (half_p - r_corner), 0)² ≤ r_corner²
                // (r_corner = 0 reduces to |u| ≤ half_m AND |v| ≤ half_p)
                double du = max(fabs(u) - inner_m, 0.0);
                double dv = max(fabs(v) - inner_p, 0.0);
                inside = du*du + dv*dv <= rc2 + eps;  // eps in Ang^2, same idea
            }
            if (inside) {
                for (int i = 0; i < 3; i++)
                    singleLoopon[i] = _compress[i]*(_u[i]*u + _v[i]*v) + _center[i];
                _loopons.push_back(singleLoopon);
            }
        }
    }
}

// =============================================================================
// snapToMidPlane
// BCC sites of generateBCC(): o + a0*Z^3 and o + a0*(Z+1/2)^3, o = -a0*(nx,ny,nz)/2.
// For coprime integers (h,k,l) their projections on n/|n| are
//   s = o.n/|n| + (a0/|n|) * (Z  U  Z + (h+k+l)/2),
// i.e. atomic planes with spacing d = a0/|n| (h+k+l even) or a0/(2|n|) (odd).
// =============================================================================
vector<double> snapToMidPlane(const vector<double>& center, const vector<double>& n,
                              double a0, int nx, int ny, int nz, double& spacing)
{
    long h[3];
    for (int i = 0; i < 3; i++) {
        h[i] = lround(n[i]);
        if (fabs(n[i] - h[i]) > 1e-9) {
            cerr << "ERROR: loop normal n must be integer Miller indices\n";
            abort();
        }
    }
    long g = 0;
    for (int i = 0; i < 3; i++) {           // gcd of |h|, |k|, |l|
        long a = labs(h[i]);
        while (a) { long t = g % a; g = a; a = t; }
    }
    if (g == 0) { cerr << "ERROR: loop normal n is zero\n"; abort(); }
    for (int i = 0; i < 3; i++) h[i] /= g;

    const double norm = sqrt(double(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]));
    const double nh[3] = { h[0]/norm, h[1]/norm, h[2]/norm };
    spacing = (labs(h[0] + h[1] + h[2]) % 2 == 0) ? a0 / norm : a0 / (2. * norm);

    const double o[3] = { -0.5*a0*nx, -0.5*a0*ny, -0.5*a0*nz };
    const double s0 = o[0]*nh[0] + o[1]*nh[1] + o[2]*nh[2];
    const double s  = center[0]*nh[0] + center[1]*nh[1] + center[2]*nh[2];
    // A centre exactly on an atomic plane goes to the mid-plane on the +n side
    const double sMid = s0 + (floor((s - s0) / spacing + 1e-9) + 0.5) * spacing;

    vector<double> out(3);
    for (int i = 0; i < 3; i++) out[i] = center[i] + (sMid - s) * nh[i];
    return out;
}

// =============================================================================
// generateSurface
// dS = area of one grid cell: da^2 (square grid), da^2 sin(angle) for the
// oblique grid of a parallelogram
// =============================================================================
void loop::generateSurface()
{
    _dS = _da * _da;
    if (_shape == LoopShape::PARALLELOGRAM) {
        const double c[3] = { _u[1]*_v[2] - _u[2]*_v[1],
                              _u[2]*_v[0] - _u[0]*_v[2],
                              _u[0]*_v[1] - _u[1]*_v[0] };
        _dS *= sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
    }
}

// =============================================================================
// Getters
// =============================================================================
LoopShape           loop::getShape()            const { return _shape; }
double              loop::getSurface()          const { return _dS; }
int                 loop::getNumberOfLoopons()  const { return (int)_loopons.size(); }
vector<double>      loop::getOneLoopon(int nb)  const { return _loopons[nb]; }
vector<double>      loop::getCenter()           const { return _center; }
vector<double>      loop::getM()                const { return _m; }
vector<double>      loop::getN()                const { return _n; }
vector<double>      loop::getP()                const { return _p; }
vector<double>      loop::getBurgers()          const { return _burgers; }

#endif
