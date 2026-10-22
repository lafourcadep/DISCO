// =============================================================================
// loop_deform.cpp — implementation of box / atom deformation from dislocation
//                   loop insertion (see loop_deform.hpp for the full API).
// =============================================================================

#include "loop_deform.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>   // std::max

// -----------------------------------------------------------------------------
// Internal 3×3 helpers
// -----------------------------------------------------------------------------

static double det3(const Mat3& M)
{
    return  M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
}

static Mat3 matmul3(const Mat3& A, const Mat3& B)
{
    Mat3 C = {};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// -----------------------------------------------------------------------------
// Jacobi eigenvalue algorithm for a 3×3 symmetric matrix A.
//
//   On output:  lam[i]   = eigenvalues  (not necessarily sorted)
//               V[*][i]  = eigenvectors (columns of V)
//
//   Invariant:  A_input  =  V · diag(lam) · V^T
// -----------------------------------------------------------------------------
static void jacobi3(Mat3 A, double lam[3], Mat3& V)
{
    // Initialise V = I
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            V[i][j] = (i == j) ? 1.0 : 0.0;

    for (int iter = 0; iter < 300; iter++) {

        // Find (p, q) with p < q that has the largest |off-diagonal| entry
        int p = 0, q = 1;
        double maxOff = std::abs(A[0][1]);
        for (int i = 0; i < 3; i++)
            for (int j = i+1; j < 3; j++)
                if (std::abs(A[i][j]) > maxOff) {
                    maxOff = std::abs(A[i][j]);
                    p = i; q = j;
                }

        if (maxOff < 1e-15) break;   // converged

        // Compute Givens angle that zeros A[p][q]
        double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t   = (tau >= 0.0)
                   ?  1.0 / (tau + std::sqrt(1.0 + tau*tau))
                   :  1.0 / (tau - std::sqrt(1.0 + tau*tau));
        double c   = 1.0 / std::sqrt(1.0 + t*t);
        double s   = t * c;

        // Update A ← G^T · A · G  (G = Givens rotation in plane (p,q))
        double App = A[p][p] - t * A[p][q];
        double Aqq = A[q][q] + t * A[p][q];
        A[p][p] = App;
        A[q][q] = Aqq;
        A[p][q] = A[q][p] = 0.0;
        for (int k = 0; k < 3; k++) {
            if (k == p || k == q) continue;
            double Akp = c * A[k][p] - s * A[k][q];
            double Akq = s * A[k][p] + c * A[k][q];
            A[k][p] = A[p][k] = Akp;
            A[k][q] = A[q][k] = Akq;
        }

        // Accumulate eigenvectors: V ← V · G
        for (int k = 0; k < 3; k++) {
            double Vkp = c * V[k][p] - s * V[k][q];
            double Vkq = s * V[k][p] + c * V[k][q];
            V[k][p] = Vkp;
            V[k][q] = Vkq;
        }
    }

    for (int i = 0; i < 3; i++) lam[i] = A[i][i];
}

// =============================================================================
// Public API
// =============================================================================

Mat3 computeU(const std::vector<LoopDesc>& loops, double V_box)
{
    // ------------------------------------------------------------------
    // 1.  Accumulate plastic distortion β[i,j] = Σ_k b_k[i]*n_k[j]*A_k/V
    // ------------------------------------------------------------------
    Mat3 beta = {};
    for (const LoopDesc& lp : loops) {
        // Normalise n
        double nx = lp.n[0], ny = lp.n[1], nz = lp.n[2];
        double nn = std::sqrt(nx*nx + ny*ny + nz*nz);
        if (nn > 0.0) { nx /= nn; ny /= nn; nz /= nn; }
        const double n3[3] = {nx, ny, nz};

        const double fac = lp.A / V_box;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                beta[i][j] += lp.b[i] * n3[j] * fac;
    }

    // ------------------------------------------------------------------
    // 2.  Deformation gradient  F = I + β
    // ------------------------------------------------------------------
    Mat3 F = {};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            F[i][j] = beta[i][j] + (i == j ? 1.0 : 0.0);

    // ------------------------------------------------------------------
    // 3.  C = F^T · F  (symmetric, positive-definite)
    // ------------------------------------------------------------------
    Mat3 C = {};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                C[i][j] += F[k][i] * F[k][j];

    // ------------------------------------------------------------------
    // 4.  Eigendecompose  C = V · diag(λ) · V^T
    // ------------------------------------------------------------------
    double lam[3];
    Mat3 V;
    jacobi3(C, lam, V);

    // ------------------------------------------------------------------
    // 5.  U = V · diag(sqrt(λ)) · V^T   (right stretch, symmetric)
    // ------------------------------------------------------------------
    Mat3 U = {};
    for (int i = 0; i < 3; i++) {
        double sigma = std::sqrt(std::max(0.0, lam[i]));
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                U[r][c] += sigma * V[r][i] * V[c][i];
    }

    // ------------------------------------------------------------------
    // 6.  Sanity check: volume conservation  (det U ≈ 1 for shear loops)
    // ------------------------------------------------------------------
    double vol_ratio = det3(U);
    if (std::abs(vol_ratio - 1.0) > 1e-4)
        std::cerr << "WARNING [computeU]: volume changed by "
                  << (vol_ratio - 1.0) * 100.0 << "% after loop insertion."
                  << "  Check that b·n = 0 for all shear loops.\n";

    return U;
}

Mat3 applyLoopDeformation(
    const std::vector<LoopDesc>& loops,
    Mat3&                        H,
    std::vector<double>&         atomX,
    std::vector<double>&         atomY,
    std::vector<double>&         atomZ,
    int                          nAtoms)
{
    const double V = std::abs(det3(H));
    const Mat3   U = computeU(loops, V);

    // H' = U · H
    H = matmul3(U, H);

    // r'_k = U · r_k
    for (int k = 0; k < nAtoms; k++) {
        const double x = atomX[k], y = atomY[k], z = atomZ[k];
        atomX[k] = U[0][0]*x + U[0][1]*y + U[0][2]*z;
        atomY[k] = U[1][0]*x + U[1][1]*y + U[1][2]*z;
        atomZ[k] = U[2][0]*x + U[2][1]*y + U[2][2]*z;
    }

    return U;
}
