// =============================================================================
// Periodic Voronoi polycrystal — see polycrystal.hpp
// =============================================================================
#include "polycrystal.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

static inline double wrapInBox(double v, double L)
{
    v -= L * floor(v / L + 0.5);              // [-L/2, L/2]
    if (v >= 0.5 * L) v -= L;                 // half-open
    return v;
}

static void rotationFromRows(Orientation& o)
{
    o.active = true;
    for (int a = 0; a < 3; a++)
        for (int k = 0; k < 3; k++) o.dir[a][k] = 0;   // not integer directions
}

Orientation randomOrientation(PortableRandom& rng)
{
    // Shoemake: uniform unit quaternion (w, x, y, z)
    const double u1 = rng.uniform(), u2 = rng.uniform(), u3 = rng.uniform();
    const double a = sqrt(1. - u1), b = sqrt(u1);
    const double x = a * sin(2. * M_PI * u2), y = a * cos(2. * M_PI * u2);
    const double z = b * sin(2. * M_PI * u3), w = b * cos(2. * M_PI * u3);
    // rotation matrix of the quaternion = crystal -> lab
    Orientation o;
    const double M[3][3] = {
        { 1 - 2*(y*y + z*z), 2*(x*y - z*w),     2*(x*z + y*w) },
        { 2*(x*y + z*w),     1 - 2*(x*x + z*z), 2*(y*z - x*w) },
        { 2*(x*z - y*w),     2*(y*z + x*w),     1 - 2*(x*x + y*y) } };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) o.R[i][j] = M[i][j];
    rotationFromRows(o);
    return o;
}

Orientation orientationFromEuler(double phi1, double Phi, double phi2)
{
    const double d = M_PI / 180.;
    const double c1 = cos(phi1 * d), s1 = sin(phi1 * d), c = cos(Phi * d), s = sin(Phi * d),
                 c2 = cos(phi2 * d), s2 = sin(phi2 * d);
    const double g[3][3] = {
        {  c1*c2 - s1*s2*c,  s1*c2 + c1*s2*c, s2*s },
        { -c1*s2 - s1*c2*c, -s1*s2 + c1*c2*c, c2*s },
        {  s1*s,            -c1*s,            c    } };
    Orientation o;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) o.R[i][j] = g[j][i];   // R = g^T
    rotationFromRows(o);
    return o;
}

void eulerFromOrientation(const Orientation& o, double& phi1, double& Phi, double& phi2)
{
    // g = R^T
    const double g33 = o.R[2][2], g31 = o.R[0][2], g32 = o.R[1][2], g13 = o.R[2][0], g23 = o.R[2][1];
    const double r = 180. / M_PI;
    Phi = acos(max(-1., min(1., g33))) * r;
    if (fabs(sin(Phi / r)) > 1e-9) {
        phi1 = atan2(g31, -g32) * r;
        phi2 = atan2(g13, g23) * r;
    } else {                                   // Phi = 0 or 180: only phi1 + phi2 defined
        phi1 = atan2(o.R[1][0], o.R[0][0]) * r;   // g12, g11
        phi2 = 0.;
    }
    if (phi1 < 0.) phi1 += 360.;
    if (phi2 < 0.) phi2 += 360.;
}

double nearestNeighbourDistance(Lattice l, double a0)
{
    return l == Lattice::FCC ? a0 / sqrt(2.) : a0 * sqrt(3.) / 2.;
}

vector<Grain> makeGrains(const PolycrystalSpec& spec)
{
    PortableRandom rng(spec.seed);
    const int n = spec.list.empty() ? spec.nGrains : (int)spec.list.size();
    vector<Grain> grains(n);
    for (int g = 0; g < n; g++) {
        const bool hasC = !spec.list.empty() && spec.hasCenter[g];
        const bool hasO = !spec.list.empty() && spec.hasOrientation[g];
        for (int k = 0; k < 3; k++)
            grains[g].seed[k] = hasC ? wrapInBox(spec.list[g].seed[k], spec.box[k])
                                     : (rng.uniform() - 0.5) * spec.box[k];
        grains[g].o = hasO ? spec.list[g].o : randomOrientation(rng);
    }
    return grains;
}

int grainOf(const vector<Grain>& grains, const double box[3], const double x[3])
{
    int best = -1;
    double bd = 1e300;
    for (size_t g = 0; g < grains.size(); g++) {
        double d2 = 0.;
        for (int k = 0; k < 3; k++) {
            double r = x[k] - grains[g].seed[k];
            r -= box[k] * floor(r / box[k] + 0.5);
            d2 += r * r;
        }
        if (d2 < bd - 1e-9) { bd = d2; best = (int)g; }
    }
    return best;
}

namespace {
struct Competitor {
    double p[3];       // seed image position (lab, unwrapped)
    double dist;       // distance to the grain's own seed
    int    key[4];     // (grain, kx, ky, kz) for tie breaking
};

bool keyLess(const int a[4], const int b[4])
{
    for (int i = 0; i < 4; i++)
        if (a[i] != b[i]) return a[i] < b[i];
    return false;
}
}

PolycrystalResult buildPolycrystal(Lattice l, double a0, const PolycrystalSpec& spec,
                                   const vector<Grain>& grainsIn)
{
    static const double bcc[2][3] = { {0, 0, 0}, {0.5, 0.5, 0.5} };
    static const double fcc[4][3] = { {0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5} };
    const double (*f)[3] = (l == Lattice::FCC) ? fcc : bcc;
    const int nb = atomsPerCell(l);
    const double* L = spec.box;
    const double tieEps = 1e-8;                // Angstrom^2

    PolycrystalResult res;
    res.grains = grainsIn;
    const int G = (int)res.grains.size();

    for (int g = 0; g < G; g++) {
        const Grain& gr = res.grains[g];
        const double* s = gr.seed;
        // Competing seed images, nearest first
        vector<Competitor> comp;
        for (int j = 0; j < G; j++)
            for (int kx = -1; kx <= 1; kx++)
                for (int ky = -1; ky <= 1; ky++)
                    for (int kz = -1; kz <= 1; kz++) {
                        if (j == g && kx == 0 && ky == 0 && kz == 0) continue;
                        Competitor c;
                        c.p[0] = res.grains[j].seed[0] + kx * L[0];
                        c.p[1] = res.grains[j].seed[1] + ky * L[1];
                        c.p[2] = res.grains[j].seed[2] + kz * L[2];
                        c.dist = sqrt((c.p[0]-s[0])*(c.p[0]-s[0]) + (c.p[1]-s[1])*(c.p[1]-s[1])
                                      + (c.p[2]-s[2])*(c.p[2]-s[2]));
                        c.key[0] = j; c.key[1] = kx; c.key[2] = ky; c.key[3] = kz;
                        comp.push_back(c);
                    }
        sort(comp.begin(), comp.end(), [](const Competitor& a, const Competitor& b) { return a.dist < b.dist; });
        const int ownKey[4] = { g, 0, 0, 0 };

        // Crystal-frame range of the lab box s +/- L/2
        double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
        for (int c = 0; c < 8; c++) {
            const double x[3] = { (c & 1 ? 0.5 : -0.5) * L[0], (c & 2 ? 0.5 : -0.5) * L[1],
                                  (c & 4 ? 0.5 : -0.5) * L[2] };
            for (int k = 0; k < 3; k++) {
                const double ck = gr.o.R[0][k]*x[0] + gr.o.R[1][k]*x[1] + gr.o.R[2][k]*x[2];
                lo[k] = min(lo[k], ck);
                hi[k] = max(hi[k], ck);
            }
        }
        const long i0 = (long)floor(lo[0] / a0) - 1, i1 = (long)ceil(hi[0] / a0) + 1;
        const long j0 = (long)floor(lo[1] / a0) - 1, j1 = (long)ceil(hi[1] / a0) + 1;
        const long k0 = (long)floor(lo[2] / a0) - 1, k1 = (long)ceil(hi[2] / a0) + 1;
        const long ni = i1 - i0 + 1;

        vector<vector<double>> buf(ni);        // per crystal plane i: kept x, y, z (deterministic merge)
        long tested = 0;
#pragma omp parallel for schedule(dynamic) reduction(+:tested)
        for (long ii = 0; ii < ni; ii++) {
            const long i = i0 + ii;
            vector<double>& out = buf[ii];
            for (long j = j0; j <= j1; j++)
                for (long k = k0; k <= k1; k++)
                    for (int b = 0; b < nb; b++) {
                        const double c[3] = { (i + f[b][0]) * a0, (j + f[b][1]) * a0, (k + f[b][2]) * a0 };
                        double r[3];
                        bool in = true;
                        for (int d = 0; d < 3 && in; d++) {
                            r[d] = gr.o.R[d][0]*c[0] + gr.o.R[d][1]*c[1] + gr.o.R[d][2]*c[2];
                            in = fabs(r[d]) <= 0.5 * L[d] + 1e-6;
                        }
                        if (!in) continue;
                        tested++;
                        const double d02 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
                        const double reach = 2. * sqrt(d02) + 1e-6;
                        const double x[3] = { s[0] + r[0], s[1] + r[1], s[2] + r[2] };
                        bool mine = true;
                        for (const Competitor& cp : comp) {
                            if (cp.dist > reach) break;   // |x - p| >= dist - |x - s| > |x - s|
                            const double e[3] = { x[0] - cp.p[0], x[1] - cp.p[1], x[2] - cp.p[2] };
                            const double d2 = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
                            if (d2 < d02 - tieEps || (fabs(d2 - d02) <= tieEps && keyLess(cp.key, ownKey))) {
                                mine = false;
                                break;
                            }
                        }
                        if (!mine) continue;
                        out.push_back(wrapInBox(x[0], L[0]));
                        out.push_back(wrapInBox(x[1], L[1]));
                        out.push_back(wrapInBox(x[2], L[2]));
                    }
        }
        res.candidates += tested;
        for (const auto& out : buf)
            for (size_t q = 0; q < out.size(); q += 3) {
                res.x.push_back(out[q]); res.y.push_back(out[q + 1]); res.z.push_back(out[q + 2]);
                res.grain.push_back(g);
            }
    }

    // ---- GB overlap removal: pairs closer than r_min lose one atom ----
    const double rmin = spec.gbMinDistance * nearestNeighbourDistance(l, a0);
    const long N = (long)res.x.size();
    int nc[3];
    for (int d = 0; d < 3; d++) nc[d] = max(1, (int)floor(L[d] / rmin));
    const bool brute = nc[0] < 3 || nc[1] < 3 || nc[2] < 3;
    auto cellOf = [&](long a, int d) {
        const double v = (d == 0 ? res.x[a] : d == 1 ? res.y[a] : res.z[a]) + 0.5 * L[d];
        int c = (int)floor(v / L[d] * nc[d]);
        return min(max(c, 0), nc[d] - 1);
    };
    vector<long> head, next(N, -1);
    if (!brute) {
        head.assign((size_t)nc[0] * nc[1] * nc[2], -1);
        for (long a = N - 1; a >= 0; a--) {
            const long c = ((long)cellOf(a, 0) * nc[1] + cellOf(a, 1)) * nc[2] + cellOf(a, 2);
            next[a] = head[c];
            head[c] = a;
        }
    }
    auto dist2 = [&](long a, long b) {
        double e[3] = { res.x[b] - res.x[a], res.y[b] - res.y[a], res.z[b] - res.z[a] };
        double s2 = 0.;
        for (int d = 0; d < 3; d++) { e[d] -= L[d] * floor(e[d] / L[d] + 0.5); s2 += e[d] * e[d]; }
        return s2;
    };
    vector<pair<long, long>> pairs;
    const double rmin2 = rmin * rmin;
#pragma omp parallel
    {
        vector<pair<long, long>> mine;
#pragma omp for schedule(dynamic, 4096)
        for (long a = 0; a < N; a++) {
            if (brute) {
                for (long b = a + 1; b < N; b++)
                    if (dist2(a, b) < rmin2) mine.push_back({ a, b });
                continue;
            }
            const int ca[3] = { cellOf(a, 0), cellOf(a, 1), cellOf(a, 2) };
            for (int dx = -1; dx <= 1; dx++)
                for (int dy = -1; dy <= 1; dy++)
                    for (int dz = -1; dz <= 1; dz++) {
                        const int cx = (ca[0] + dx + nc[0]) % nc[0], cy = (ca[1] + dy + nc[1]) % nc[1],
                                  cz = (ca[2] + dz + nc[2]) % nc[2];
                        for (long b = head[((long)cx * nc[1] + cy) * nc[2] + cz]; b >= 0; b = next[b])
                            if (b > a && dist2(a, b) < rmin2) mine.push_back({ a, b });
                    }
        }
#pragma omp critical
        pairs.insert(pairs.end(), mine.begin(), mine.end());
    }
    sort(pairs.begin(), pairs.end());
    pairs.erase(unique(pairs.begin(), pairs.end()), pairs.end());   // nc < 3 never reached: brute
    vector<char> removed(N, 0);
    for (const auto& pr : pairs) {
        const long a = pr.first, b = pr.second;
        if (removed[a] || removed[b]) continue;
        const long drop = (res.grain[b] > res.grain[a] || (res.grain[b] == res.grain[a])) ? b : a;
        removed[drop] = 1;
        res.removedGB++;
    }
    long k = 0;
    for (long a = 0; a < N; a++) {
        if (removed[a]) continue;
        res.x[k] = res.x[a]; res.y[k] = res.y[a]; res.z[k] = res.z[a]; res.grain[k] = res.grain[a];
        k++;
    }
    res.x.resize(k); res.y.resize(k); res.z.resize(k); res.grain.resize(k);
    for (auto& gr : res.grains) gr.atoms = 0;
    for (long a = 0; a < k; a++) res.grains[res.grain[a]].atoms++;
    return res;
}

// =============================================================================
// Cells as half-spaces, loop clipping
// =============================================================================
vector<HalfSpace> cellHalfSpaces(const vector<Grain>& grains, const double box[3], int k)
{
    vector<HalfSpace> hs;
    const double* s = grains[k].seed;
    for (size_t j = 0; j < grains.size(); j++)
        for (int kx = -1; kx <= 1; kx++)
            for (int ky = -1; ky <= 1; ky++)
                for (int kz = -1; kz <= 1; kz++) {
                    if ((int)j == k && kx == 0 && ky == 0 && kz == 0) continue;
                    HalfSpace h;
                    const double p[3] = { grains[j].seed[0] + kx * box[0], grains[j].seed[1] + ky * box[1],
                                          grains[j].seed[2] + kz * box[2] };
                    double a2 = 0., as = 0.;
                    for (int d = 0; d < 3; d++) {
                        h.a[d] = p[d] - s[d];
                        a2 += h.a[d] * h.a[d];
                        as += h.a[d] * s[d];
                    }
                    h.c = as + 0.5 * a2;
                    hs.push_back(h);
                }
    return hs;
}

bool insideCell(const vector<HalfSpace>& hs, const double x[3], double tol)
{
    for (const HalfSpace& h : hs)
        if (h.a[0]*x[0] + h.a[1]*x[1] + h.a[2]*x[2] > h.c + tol) return false;
    return true;
}

double polygonArea(const vector<vector<double>>& P)
{
    double c[3] = { 0., 0., 0. };
    const size_t N = P.size();
    for (size_t i = 0; i < N; i++) {
        const vector<double>& p = P[i];
        const vector<double>& q = P[(i + 1) % N];
        c[0] += p[1]*q[2] - p[2]*q[1];
        c[1] += p[2]*q[0] - p[0]*q[2];
        c[2] += p[0]*q[1] - p[1]*q[0];
    }
    return 0.5 * sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
}

vector<vector<double>> clipPolygon(const vector<vector<double>>& poly, const vector<HalfSpace>& hs)
{
    vector<vector<double>> P = poly;
    for (const HalfSpace& h : hs) {
        if (P.empty()) break;
        const double an = sqrt(h.a[0]*h.a[0] + h.a[1]*h.a[1] + h.a[2]*h.a[2]);
        auto f = [&](const vector<double>& x) {    // signed distance, > 0 outside
            return (h.a[0]*x[0] + h.a[1]*x[1] + h.a[2]*x[2] - h.c) / an;
        };
        bool anyOut = false;
        for (const auto& x : P) if (f(x) > 1e-12) { anyOut = true; break; }
        if (!anyOut) continue;
        vector<vector<double>> Q;
        const size_t N = P.size();
        for (size_t i = 0; i < N; i++) {
            const vector<double>& A = P[i];
            const vector<double>& B = P[(i + 1) % N];
            const double fa = f(A), fb = f(B);
            const bool ina = fa <= 1e-12, inb = fb <= 1e-12;
            if (ina) Q.push_back(A);
            if (ina != inb) {
                const double t = fa / (fa - fb);
                Q.push_back({ A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1]), A[2] + t * (B[2] - A[2]) });
            }
        }
        P.swap(Q);
    }
    // Clean up: merge near-duplicates, drop collinear vertices
    bool changed = true;
    while (changed && P.size() >= 3) {
        changed = false;
        const size_t N = P.size();
        for (size_t i = 0; i < N && !changed; i++) {
            const vector<double>& A = P[(i + N - 1) % N];
            const vector<double>& B = P[i];
            const vector<double>& C = P[(i + 1) % N];
            const double e1[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] }, e2[3] = { C[0]-B[0], C[1]-B[1], C[2]-B[2] };
            const double l1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
            const double l2 = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
            const double cr[3] = { e1[1]*e2[2] - e1[2]*e2[1], e1[2]*e2[0] - e1[0]*e2[2], e1[0]*e2[1] - e1[1]*e2[0] };
            const double sinT = (l1 > 0. && l2 > 0.) ? sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]) / (l1 * l2) : 0.;
            if (l1 < 1e-6 || sinT < 1e-9) {        // B duplicates A, or A-B-C collinear
                P.erase(P.begin() + i);
                changed = true;
            }
        }
    }
    if (P.size() < 3) P.clear();
    return P;
}
