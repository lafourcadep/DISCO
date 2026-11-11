// =============================================================================
// Coherent twin lamellae — see twin.hpp
// =============================================================================
#include "twin.hpp"

#include <cmath>
#include <iostream>

using namespace std;

static void twinFail(const string& msg)
{
    cerr << "ERROR in twins: " << msg << "\n";
    abort();
}

static inline double dot3(const double a[3], const double b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// Unwrap x around the seed (polycrystal) or leave it (single crystal)
static void unwrap(const double x[3], const double* seed, const double box[3], double out[3])
{
    for (int k = 0; k < 3; k++) {
        if (!seed) { out[k] = x[k]; continue; }
        double r = x[k] - seed[k];
        r -= box[k] * floor(r / box[k] + 0.5);
        out[k] = seed[k] + r;
    }
}

Twin makeTwin(const TwinSpec& spec, Lattice l, double a0, const Orientation& parent,
              const double origin[3], const double* seed, const double box[3], double& snapShift)
{
    Twin t;
    t.grain = spec.grain;
    t.parent = parent;
    if (!isTwinPlane(l, spec.n))
        twinFail(string("n = [") + to_string(spec.n[0]) + " " + to_string(spec.n[1]) + " " + to_string(spec.n[2])
                 + "] is not a twin plane (" + (l == Lattice::FCC ? "FCC: {111}" : "BCC: {112}") + ")");
    long g = 0;
    for (int i = 0; i < 3; i++) { t.h[i] = lround(spec.n[i]); long a = labs(t.h[i]); while (a) { long r = g % a; g = a; a = r; } }
    for (int i = 0; i < 3; i++) t.h[i] /= g;
    const double hn = sqrt(double(t.h[0]*t.h[0] + t.h[1]*t.h[1] + t.h[2]*t.h[2]));
    const vector<double> nC = { t.h[0] / hn, t.h[1] / hn, t.h[2] / hn };
    const vector<double> nL = rotateToLab(parent, nC);
    for (int k = 0; k < 3; k++) t.nLab[k] = nL[k];
    t.d = latticePlaneSpacing(l, t.h, a0);
    const double period = latticePeriod(l, t.h, a0);
    const long q = lround(period / t.d);                  // layers per period along n
    const long m0 = (q % 2 == 0) ? q / 2 : q;             // smallest m with 2 m = 0 mod q
    t.m = (int)max(m0, m0 * lround(spec.thickness / (t.d * m0)));

    // A lattice vector v with v.n = d (crystal: a0/2 (h,k,l), parity of the lattice)
    bool found = false;
    for (int a = -2; a <= 2 && !found; a++)
        for (int b = -2; b <= 2 && !found; b++)
            for (int c = -2; c <= 2 && !found; c++) {
                if (!isLatticeTranslation(l, { double(a), double(b), double(c) })) continue;
                const double v[3] = { 0.5 * a0 * a, 0.5 * a0 * b, 0.5 * a0 * c };
                if (fabs(v[0]*nC[0] + v[1]*nC[1] + v[2]*nC[2] - t.d) > 1e-9) continue;
                const vector<double> vl = rotateToLab(parent, { v[0], v[1], v[2] });
                for (int k = 0; k < 3; k++) t.vLab[k] = vl[k];
                found = true;
            }
    if (!found) twinFail("no lattice vector with one plane spacing along n");

    // First TB: the atomic plane nearest to the centre (unwrapped around the seed)
    double c[3];
    unwrap(spec.center.data(), seed, box, c);
    const double s0 = dot3(origin, t.nLab), s = dot3(c, t.nLab);
    const long K = lround((s - s0) / t.d);
    snapShift = s0 + K * t.d - s;
    // O: the parent site of that plane nearest to the snapped centre (kept close to the grain,
    // so points built from O unwrap consistently around the seed)
    {
        double cp[3];
        for (int k = 0; k < 3; k++) cp[k] = c[k] + snapShift * t.nLab[k];
        const vector<double> cc = rotateToCrystal(parent, { cp[0] - origin[0], cp[1] - origin[1], cp[2] - origin[2] });
        const double cq[3] = { cc[0], cc[1], cc[2] };
        double cz[3];
        nearestLatticeSite(l, cq, a0, cz);
        const vector<double> zl = rotateToLab(parent, { cz[0], cz[1], cz[2] });
        for (int k = 0; k < 3; k++) t.O[k] = origin[k] + zl[k];
        const long off = lround((dot3(t.O, t.nLab) - (s0 + K * t.d)) / t.d);   // nearest site may be a plane off
        for (int k = 0; k < 3; k++) t.O[k] -= off * t.vLab[k];
        if (fabs(dot3(t.O, t.nLab) - (s0 + K * t.d)) > 1e-6) twinFail("internal: TB origin not on its plane");
    }

    // Single crystal: n along a box axis, box a whole number of planes along it
    if (!seed) {
        for (int k = 0; k < 3; k++) if (fabs(fabs(t.nLab[k]) - 1.) < 1e-9) t.axis = k;
        if (t.axis < 0)
            twinFail("in a single crystal the twin plane must be normal to a box axis (use Orientation, "
                     "e.g. Z = the twin normal)");
        const double nl = box[t.axis] / t.d;
        t.layersPerBox = (int)lround(nl);
        if (fabs(nl - t.layersPerBox) > 1e-6 || t.layersPerBox % q != 0)
            twinFail("the box is not a whole number of lattice periods along the twin normal");
        if (t.m >= t.layersPerBox)
            twinFail("the lamella is as thick as the box: no matrix left");
    }

    // Twin frame: R_twin = R_parent * halfTurn(n_crystal)
    double H[3][3];
    const double u[3] = { nC[0], nC[1], nC[2] };
    halfTurn(u, H);
    t.twin = parent;
    t.twin.active = true;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double r = 0.;
            for (int k = 0; k < 3; k++) r += parent.R[i][k] * H[k][j];
            t.twin.R[i][j] = r;
        }
    for (int a = 0; a < 3; a++) for (int k = 0; k < 3; k++) t.twin.dir[a][k] = 0;

    // Per-layer in-plane shifts t_j
    auto mirror = [&](const double x[3], double out[3]) {
        double e[3] = { x[0] - t.O[0], x[1] - t.O[1], x[2] - t.O[2] };
        const double sn = dot3(e, t.nLab);
        for (int k = 0; k < 3; k++) out[k] = x[k] - 2. * sn * t.nLab[k];
    };
    t.shift.assign(t.m + 1, vector<double>(3, 0.));
    for (int j = 1; j < t.m; j++) {
        double xr[3], q3[3], z[3], xp[3];
        for (int k = 0; k < 3; k++) xr[k] = t.O[k] + j * t.vLab[k];
        mirror(xr, q3);
        const vector<double> cq = rotateToCrystal(parent, { q3[0] - origin[0], q3[1] - origin[1], q3[2] - origin[2] });
        const double cc[3] = { cq[0], cq[1], cq[2] };
        double cz[3];
        nearestLatticeSite(l, cc, a0, cz);
        const vector<double> zl = rotateToLab(parent, { cz[0], cz[1], cz[2] });
        for (int k = 0; k < 3; k++) z[k] = origin[k] + zl[k];
        mirror(z, xp);
        for (int k = 0; k < 3; k++) t.shift[j][k] = xp[k] - xr[k];
        const double along = dot3(t.shift[j].data(), t.nLab);
        if (fabs(along) > 1e-6) twinFail("internal: twin shift not in the plane");
        // A parent atom sits at a hollow of the mirrored net: several in-plane lattice
        // equivalents of t_j are equally short. Pick the shortest, ties broken by the
        // lab components, so the result does not depend on rounding or on the choice of O.
        vector<double> best = t.shift[j];
        double bestN = dot3(best.data(), best.data());
        for (int a = -2; a <= 2; a++)
            for (int b = -2; b <= 2; b++)
                for (int c = -2; c <= 2; c++) {
                    if (!isLatticeTranslation(l, { double(a), double(b), double(c) })) continue;
                    const vector<double> w = rotateToLab(parent, { 0.5 * a0 * a, 0.5 * a0 * b, 0.5 * a0 * c });
                    const double wl[3] = { w[0], w[1], w[2] };
                    if (fabs(dot3(wl, t.nLab)) > 1e-9) continue;          // in-plane translations only
                    const vector<double> cand = { t.shift[j][0] + w[0], t.shift[j][1] + w[1], t.shift[j][2] + w[2] };
                    const double cn = dot3(cand.data(), cand.data());
                    bool better = cn < bestN - 1e-9;
                    if (!better && fabs(cn - bestN) <= 1e-9)
                        for (int q = 0; q < 3; q++) {
                            if (cand[q] < best[q] - 1e-9) { better = true; break; }
                            if (cand[q] > best[q] + 1e-9) break;
                        }
                    if (better) { best = cand; bestN = cn; }
                }
        t.shift[j] = best;
    }
    return t;
}

long twinLayer(const Twin& t, const double x[3], const double* seed, const double box[3])
{
    double xu[3];
    unwrap(x, seed, box, xu);
    const double e[3] = { xu[0] - t.O[0], xu[1] - t.O[1], xu[2] - t.O[2] };
    long j = lround(dot3(e, t.nLab) / t.d);
    if (t.axis >= 0) {                         // periodic along the box axis
        j %= t.layersPerBox;
        if (j < 0) j += t.layersPerBox;
    }
    return j;
}

void applyTwin(Twin& t, int index, vector<vector<double>>& atoms, const vector<int>& atomGrain,
               vector<int>& atomTwin, const double* seed, const double box[3])
{
    long moved = 0;
    const long N = (long)atoms.size();
#pragma omp parallel for reduction(+:moved)
    for (long i = 0; i < N; i++) {
        if (t.grain >= 0 && atomGrain[i] != t.grain) continue;
        const double x[3] = { atoms[i][0], atoms[i][1], atoms[i][2] };
        const long j = twinLayer(t, x, seed, box);
        if (j < 1 || j >= t.m) continue;
        for (int k = 0; k < 3; k++) {
            double v = atoms[i][k] + t.shift[j][k];
            v -= box[k] * floor(v / box[k] + 0.5);
            if (v >= 0.5 * box[k]) v -= box[k];
            atoms[i][k] = v;
        }
        atomTwin[i] = index + 1;
        moved++;
    }
    t.atomsMoved = moved;
}

int twinRegion(const vector<Twin>& twins, int grain, const double x[3], const double* seed, const double box[3])
{
    for (size_t i = 0; i < twins.size(); i++) {
        if (twins[i].grain != grain) continue;
        double xu[3];
        unwrap(x, seed, box, xu);
        const double e[3] = { xu[0] - twins[i].O[0], xu[1] - twins[i].O[1], xu[2] - twins[i].O[2] };
        double s = dot3(e, twins[i].nLab) / twins[i].d;
        if (twins[i].axis >= 0) {
            s = fmod(s, (double)twins[i].layersPerBox);
            if (s < 0.) s += twins[i].layersPerBox;
        }
        if (s > 0. && s < twins[i].m) return (int)i + 1;
    }
    return 0;
}
