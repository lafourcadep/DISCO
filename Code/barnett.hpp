#pragma once
// =============================================================================
// Line-integral (closed-form) displacement field of a dislocation loop made of
// straight segments, isotropic elasticity — the exact counterpart of the Kroupa
// loopon sum (kroupa.hpp), with the same nu and the same sign convention
// (displacement jump u+ - u- = -b across the cut, + on the side n points to).
//
// Burgers / de Wit form of the Volterra field of a loop C bounding a flat cut S:
//   u(x) = S1 b Omega(x) / (4 pi)
//        + sum over segments of closed-form line integrals,
// Omega = solid angle of S seen from x (fan of triangles from the loop centre,
// Van Oosterom & Strackee 1983). For a segment A -> B with unit tangent t and
// w = b x t, writing R(s) = d + s t (R = x' - x, d the part of A - x normal to t,
// s = R.t from sA to sB, R = |R|):
//   J0 = int ds / R          = ln[(R_B + s_B) / (R_A + s_A)]
//   K0 = int ds / R^3        = [s / (d^2 R)]
//   K1 = int s ds / R^3      = [-1 / R]
//   line term = w J0 (c3 - 1 / (4 pi)) - c3 (d.w) (d K0 + t K1),
//   c3 = 1 / (8 pi (1 - nu)),
// and the solid-angle term is + b Omega / (4 pi), with C running counter-clockwise
// about n (Omega > 0 on the side n points to).
// J0 and K0 use cancellation-free forms (see below). The signs of the three terms
// and the orientation of C were fixed against a direct quadrature of the Kroupa
// kernel over S: of the 16 sign / orientation combinations, only this one (and
// its overall mirror) agrees, to 2e-10 Ang (tests/test_barnett.cpp checks it on
// glide, prismatic and mixed loops of several shapes).
//
// The field is singular on C itself; DISCO's cut lies on a mid-plane between
// atomic planes, so no atom is closer to C than half a plane spacing.
// =============================================================================

#include <cmath>

#include "kroupa.hpp"   // KROUPA_HD, KROUPA_NU

KROUPA_HD double barnettLog(double x)  { return log(x); }
KROUPA_HD float  barnettLog(float x)   { return logf(x); }
KROUPA_HD double barnettAtan2(double y, double x) { return atan2(y, x); }
KROUPA_HD float  barnettAtan2(float y, float x)   { return atan2f(y, x); }

// Adds to (ux, uy, uz) the field at point p of the closed polyline
// V[0..nv-1] (V[nv-1] -> V[0] closes it), cut fanned from centre c, Burgers
// vector b. Coordinates are any consistent frame (p, V, c in the same one).
template <typename R>
KROUPA_HD void barnettLoop(R px, R py, R pz,
                           const R* vx, const R* vy, const R* vz, int nv,
                           R cx, R cy, R cz, R bx, R by, R bz,
                           R& ux, R& uy, R& uz)
{
    const R PI  = R(3.14159265358979323846);
    const R c3  = R(1.0 / (8.0 * 3.14159265358979323846 * (1.0 - KROUPA_NU)));
    const R cJ  = c3 - R(1) / (R(4) * PI);
    const R cK  = c3;

    // r1 = c - p (fan apex), shared by all triangles
    const R r1x = cx - px, r1y = cy - py, r1z = cz - pz;
    const R l1  = sqrt(r1x*r1x + r1y*r1y + r1z*r1z);

    R omega = R(0);
    R lx = R(0), ly = R(0), lz = R(0);

    R ax = vx[nv-1] - px, ay = vy[nv-1] - py, az = vz[nv-1] - pz;   // A - p
    R la = sqrt(ax*ax + ay*ay + az*az);
    for (int k = 0; k < nv; k++) {
        const R bxp = vx[k] - px, byp = vy[k] - py, bzp = vz[k] - pz;   // B - p
        const R lb  = sqrt(bxp*bxp + byp*byp + bzp*bzp);

        // ---- solid angle of triangle (c, A, B) ----
        const R num = r1x*(ay*bzp - az*byp) + r1y*(az*bxp - ax*bzp) + r1z*(ax*byp - ay*bxp);
        const R den = l1*la*lb + (r1x*ax + r1y*ay + r1z*az)*lb
                    + (r1x*bxp + r1y*byp + r1z*bzp)*la + (ax*bxp + ay*byp + az*bzp)*l1;
        omega += R(2) * barnettAtan2(num, den);

        // ---- segment A -> B ----
        R tx = bxp - ax, ty = byp - ay, tz = bzp - az;
        const R L = sqrt(tx*tx + ty*ty + tz*tz);
        if (L > R(0)) {
            tx /= L; ty /= L; tz /= L;
            const R sA = ax*tx + ay*ty + az*tz;            // R(s) = A - p + ... , s = R.t
            const R sB = sA + L;
            const R dx = ax - sA*tx, dy = ay - sA*ty, dz = az - sA*tz;
            const R d2 = dx*dx + dy*dy + dz*dz;
            const R RA = la, RB = lb;
            // J0 = ln f(B) - ln f(A), f = R + s, or d^2 / (R - s) when s < 0
            const R fA = sA >= R(0) ? RA + sA : d2 / (RA - sA);
            const R fB = sB >= R(0) ? RB + sB : d2 / (RB - sB);
            const R J0 = barnettLog(fB / fA);
            // K0 = (sB/RB - sA/RA) / d^2 without cancellation
            R K0;
            if (sA >= R(0))      K0 = R(1) / (RA*(RA + sA)) - R(1) / (RB*(RB + sB));
            else if (sB <= R(0)) K0 = R(1) / (RB*(RB - sB)) - R(1) / (RA*(RA - sA));
            else                 K0 = (sB/RB - sA/RA) / d2;
            const R K1 = R(1) / RA - R(1) / RB;
            // w = b x t
            const R wx = by*tz - bz*ty, wy = bz*tx - bx*tz, wz = bx*ty - by*tx;
            const R dw = dx*wx + dy*wy + dz*wz;
            lx += wx*J0*cJ - cK*dw*(dx*K0 + tx*K1);
            ly += wy*J0*cJ - cK*dw*(dy*K0 + ty*K1);
            lz += wz*J0*cJ - cK*dw*(dz*K0 + tz*K1);
        }
        ax = bxp; ay = byp; az = bzp; la = lb;
    }
    const R fo = omega / (R(4) * PI);
    ux += bx*fo + lx;
    uy += by*fo + ly;
    uz += bz*fo + lz;
}
