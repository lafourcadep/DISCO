#ifndef LOOP_CPP
#define LOOP_CPP

#include "loop.hpp"

#include <algorithm>

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
// In-plane loopon-grid axes when no direction is given: u = m for a glide loop
// (b in the plane). For a loop with b.n != 0, u is the projection of m on the
// plane or, for a pure prismatic loop (m // n), n x e_k with e_k the Cartesian
// axis of smallest |n_k|. Always v = u x n.
// =============================================================================
static void defaultAxes(const vector<double>& m, const vector<double>& n,
                        vector<double>& u, vector<double>& v)
{
    const double mn = m[0]*n[0] + m[1]*n[1] + m[2]*n[2];
    u = { m[0] - mn*n[0], m[1] - mn*n[1], m[2] - mn*n[2] };
    double l = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    if (l < 1e-6) {
        int k = 0;
        for (int i = 1; i < 3; i++) if (fabs(n[i]) < fabs(n[k])) k = i;
        double e[3] = { 0., 0., 0. };
        e[k] = 1.;
        u = { n[1]*e[2] - n[2]*e[1], n[2]*e[0] - n[0]*e[2], n[0]*e[1] - n[1]*e[0] };
        l = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    }
    for (double& x : u) x /= l;
    v = { u[1]*n[2] - u[2]*n[1], u[2]*n[0] - u[0]*n[2], u[0]*n[1] - u[1]*n[0] };
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
    defaultAxes(_m, _n, _u, _v);
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
        defaultAxes(_m, _n, _u, _v);
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
// Convex polygon constructor
// The vertices are projected on the orthonormal grid basis (_u, _v), put in
// counter-clockwise order in (u, v) and checked for strict convexity. With
// corner_r > 0, every side moves inward by corner_r; the inner polygon Q has the
// vertex Q_i = P_i + r (n1 + n2) / (1 + n1.n2) (n1, n2: inward unit normals of the
// two sides at P_i) and the rounded contour is { x : dist(x, Q) <= r }.
// A side of length L between interior angles a_i, a_j shrinks by
// r (cot(a_i/2) + cot(a_j/2)), hence r_max = min L / (cot(a_i/2) + cot(a_j/2)).
// =============================================================================
loop::loop(double da, const vector< vector<double> >& vertices, double corner_r,
           vector<double> center, vector<double> m, vector<double> n,
           double latticeParameter, vector<double> compress,
           vector<double> axis)
  : _shape(LoopShape::POLYGON),
    _radius(0.),
    _half_m(0.), _half_p(0.), _r_corner(0.),
    _center(center), _m(m), _n(n),
    _burgers(3, 0.), _p(3, 0.),
    _compress(compress), _a0(latticeParameter), _da(da)
{
    initCommon(_m, m, _n, n, _p, _burgers, _compress, _a0);
    if (axis.empty()) {
        defaultAxes(_m, _n, _u, _v);
    } else {
        const double an = axis[0]*_n[0] + axis[1]*_n[1] + axis[2]*_n[2];
        _u = { axis[0] - an*_n[0], axis[1] - an*_n[1], axis[2] - an*_n[2] };
        const double l = sqrt(_u[0]*_u[0] + _u[1]*_u[1] + _u[2]*_u[2]);
        for (double& x : _u) x /= l;
        _v = { _u[1]*_n[2] - _u[2]*_n[1],
               _u[2]*_n[0] - _u[0]*_n[2],
               _u[0]*_n[1] - _u[1]*_n[0] };
    }

    const int N = (int)vertices.size();
    if (N < 3) { cerr << "ERROR: a polygon loop needs at least 3 vertices\n"; abort(); }
    for (const auto& x : vertices) {
        _pu.push_back(x[0]*_u[0] + x[1]*_u[1] + x[2]*_u[2]);
        _pv.push_back(x[0]*_v[0] + x[1]*_v[1] + x[2]*_v[2]);
    }
    double area2 = 0.;
    for (int i = 0; i < N; i++) {
        const int j = (i + 1) % N;
        area2 += _pu[i]*_pv[j] - _pu[j]*_pv[i];
    }
    if (fabs(area2) < 1e-12) { cerr << "ERROR: polygon loop has zero area\n"; abort(); }
    if (area2 < 0.) {
        reverse(_pu.begin(), _pu.end());
        reverse(_pv.begin(), _pv.end());
    }

    // Unit side directions d_i = (P_{i+1} - P_i) / L_i
    vector<double> du(N), dv(N), len(N);
    for (int i = 0; i < N; i++) {
        const int j = (i + 1) % N;
        const double eu = _pu[j] - _pu[i], ev = _pv[j] - _pv[i];
        len[i] = sqrt(eu*eu + ev*ev);
        if (len[i] < 1e-9) { cerr << "ERROR: polygon loop has two equal consecutive vertices\n"; abort(); }
        du[i] = eu / len[i];
        dv[i] = ev / len[i];
    }
    // Strict convexity: every turn to the left, total turn 2 pi (rules out stars)
    vector<double> cotHalf(N);   // cot(interior angle / 2) at vertex i
    double turn = 0.;
    for (int i = 0; i < N; i++) {
        const int k = (i + N - 1) % N;
        const double sn = du[k]*dv[i] - dv[k]*du[i];
        const double cs = du[k]*du[i] + dv[k]*dv[i];
        if (sn < 1e-9) {
            cerr << "ERROR: polygon loop must be strictly convex, its vertices given in order"
                    " around the contour\n";
            abort();
        }
        turn += atan2(sn, cs);
        cotHalf[i] = (1. - cs) / sn;          // interior angle a = pi - turn angle
    }
    if (fabs(turn - 2.*M_PI) > 1e-6) {
        cerr << "ERROR: polygon loop is self-intersecting\n";
        abort();
    }

    double rMax = 1e300;
    for (int i = 0; i < N; i++)
        rMax = min(rMax, len[i] / (cotHalf[i] + cotHalf[(i + 1) % N]));
    if (corner_r < 0.) { cerr << "ERROR: Corner_radius must be >= 0\n"; abort(); }
    if (corner_r > rMax) {
        cerr << "WARNING: corner_radius " << corner_r
             << " clamped to the polygon maximum " << rMax << "\n";
        corner_r = rMax;
    }
    _r_corner = corner_r;

    _nu.resize(N); _nv.resize(N);
    for (int i = 0; i < N; i++) { _nu[i] = -dv[i]; _nv[i] = du[i]; }   // inward (CCW)
    _qu.resize(N); _qv.resize(N);
    for (int i = 0; i < N; i++) {
        const int k = (i + N - 1) % N;
        const double n1u = -dv[k], n1v = du[k], n2u = -dv[i], n2v = du[i];   // inward (CCW)
        const double f = _r_corner / (1. + n1u*n2u + n1v*n2v);
        _qu[i] = _pu[i] + f * (n1u + n2u);
        _qv[i] = _pv[i] + f * (n1v + n2v);
    }

    for (int i = 0; i < N; i++) {
        _half_m = max(_half_m, fabs(_pu[i]));
        _half_p = max(_half_p, fabs(_pv[i]));
    }
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

    const double ext_m = (_shape == LoopShape::CIRCULAR) ? _radius : _half_m;
    const double ext_p = (_shape == LoopShape::CIRCULAR) ? _radius : _half_p;
    const int steps_m = (int)ceil(ext_m / _da) + 1;
    const int steps_p = (int)ceil(ext_p / _da) + 1;

    for (int xi = -steps_m; xi <= steps_m; xi++) {
        for (int yi = -steps_p; yi <= steps_p; yi++) {
            const double u = xi * _da;
            const double v = yi * _da;
            if (insideUV(u, v)) {
                for (int i = 0; i < 3; i++)
                    singleLoopon[i] = _compress[i]*(_u[i]*u + _v[i]*v) + _center[i];
                _loopons.push_back(singleLoopon);
            }
        }
    }
}

// =============================================================================
// insideUV
// Grid points exactly on the contour: excluded for circles, included for
// rectangles, parallelograms and polygons. The tolerance makes this robust to rounding
// once sizes are given in Angstrom.
// =============================================================================
bool loop::insideUV(double u, double v) const
{
    const double eps = 1e-9;
    if (_shape == LoopShape::CIRCULAR)
        return sqrt(u*u + v*v) < _radius - eps;
    if (_shape == LoopShape::POLYGON) {
        // Inside the inner polygon Q: at least r_corner inside every side of the
        // polygon (boundary included). Uses the sides of P, so it stays valid when Q
        // degenerates to a segment or a point at the largest corner radius ...
        const int N = (int)_qu.size();
        bool inQ = true;
        for (int i = 0; i < N && inQ; i++)
            if (_nu[i]*(u - _pu[i]) + _nv[i]*(v - _pv[i]) < _r_corner - eps) inQ = false;
        if (inQ) return true;
        if (_r_corner <= 0.) return false;
        // ... or within r_corner of its boundary
        double d2 = 1e300;
        for (int i = 0; i < N; i++) {
            const int j = (i + 1) % N;
            const double eu = _qu[j] - _qu[i], ev = _qv[j] - _qv[i];
            const double l2 = eu*eu + ev*ev;
            double t = l2 > 0. ? ((u - _qu[i])*eu + (v - _qv[i])*ev) / l2 : 0.;
            t = min(max(t, 0.), 1.);
            const double wu = u - _qu[i] - t*eu, wv = v - _qv[i] - t*ev;
            d2 = min(d2, wu*wu + wv*wv);
        }
        return d2 <= _r_corner*_r_corner + eps;
    }
    // Rounded rectangle:
    //   max(|u| - (half_m - r_corner), 0)² + max(|v| - (half_p - r_corner), 0)² ≤ r_corner²
    // (r_corner = 0 reduces to |u| ≤ half_m AND |v| ≤ half_p; a parallelogram has
    //  r_corner = 0 and u, v are its oblique coordinates)
    const double du = max(fabs(u) - max(_half_m - _r_corner, 0.0), 0.0);
    const double dv = max(fabs(v) - max(_half_p - _r_corner, 0.0), 0.0);
    return du*du + dv*dv <= _r_corner*_r_corner + eps;  // eps in Ang^2, same idea
}

// =============================================================================
// contains
// Coordinates of rel in the (possibly oblique) basis _u, _v of the loop plane:
// solve the 2x2 Gram system of rel.u, rel.v (identity for orthogonal axes).
// =============================================================================
bool loop::contains(const double rel[3]) const
{
    const double ru  = rel[0]*_u[0] + rel[1]*_u[1] + rel[2]*_u[2];
    const double rv  = rel[0]*_v[0] + rel[1]*_v[1] + rel[2]*_v[2];
    const double c   = _u[0]*_v[0] + _u[1]*_v[1] + _u[2]*_v[2];   // cos(u, v)
    const double det = 1.0 - c*c;
    return insideUV((ru - c*rv) / det, (rv - c*ru) / det);
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
double              loop::getCornerRadius()     const { return _r_corner; }
int                 loop::getNumberOfVertices() const { return (int)_pu.size(); }
double              loop::getSurface()          const { return _dS; }
int                 loop::getNumberOfLoopons()  const { return (int)_loopons.size(); }
vector<double>      loop::getOneLoopon(int nb)  const { return _loopons[nb]; }
vector<double>      loop::getCenter()           const { return _center; }
vector<double>      loop::getM()                const { return _m; }
vector<double>      loop::getN()                const { return _n; }
vector<double>      loop::getP()                const { return _p; }
vector<double>      loop::getBurgers()          const { return _burgers; }

// =============================================================================
// getContourArea
// Exact area enclosed by the contour (the loopon count times dS approaches it as
// Loopon_size -> 0). Rounded polygon: area(Q) + perimeter(Q) r + pi r^2.
// =============================================================================
double loop::getContourArea() const
{
    const double r = _r_corner;
    switch (_shape) {
        case LoopShape::CIRCULAR:
            return M_PI * _radius * _radius;
        case LoopShape::RECTANGULAR:
            return 4. * _half_m * _half_p - (4. - M_PI) * r * r;
        case LoopShape::PARALLELOGRAM: {
            const double c = _u[0]*_v[0] + _u[1]*_v[1] + _u[2]*_v[2];
            return 4. * _half_m * _half_p * sqrt(max(1. - c*c, 0.));
        }
        case LoopShape::POLYGON: {
            const int N = (int)_qu.size();
            double a2 = 0., per = 0.;
            for (int i = 0; i < N; i++) {
                const int j = (i + 1) % N;
                a2  += _qu[i]*_qv[j] - _qu[j]*_qv[i];
                per += sqrt((_qu[j]-_qu[i])*(_qu[j]-_qu[i]) + (_qv[j]-_qv[i])*(_qv[j]-_qv[i]));
            }
            return 0.5 * a2 + per * r + M_PI * r * r;
        }
    }
    return 0.;
}

// =============================================================================
// regularPolygon
// =============================================================================
vector< vector<double> > regularPolygon(int sides, double R, vector<double> vertexDir,
                                        const vector<double>& m, const vector<double>& n)
{
    if (sides < 3) { cerr << "ERROR: a polygon loop needs Sides >= 3\n"; abort(); }
    const double lm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
    const double ln = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    const vector<double> mh = { m[0]/lm, m[1]/lm, m[2]/lm }, nh = { n[0]/ln, n[1]/ln, n[2]/ln };
    vector<double> u, v;
    if (vertexDir.empty()) {
        defaultAxes(mh, nh, u, v);
    } else {
        const double an = vertexDir[0]*nh[0] + vertexDir[1]*nh[1] + vertexDir[2]*nh[2];
        u = { vertexDir[0] - an*nh[0], vertexDir[1] - an*nh[1], vertexDir[2] - an*nh[2] };
        const double l = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
        for (double& x : u) x /= l;
        v = { u[1]*nh[2] - u[2]*nh[1], u[2]*nh[0] - u[0]*nh[2], u[0]*nh[1] - u[1]*nh[0] };
    }
    vector< vector<double> > out(sides, vector<double>(3));
    for (int k = 0; k < sides; k++) {
        const double th = 2. * M_PI * k / sides;
        for (int i = 0; i < 3; i++) out[k][i] = R * (cos(th) * u[i] + sin(th) * v[i]);
    }
    return out;
}

#endif
