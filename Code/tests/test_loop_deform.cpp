// =============================================================================
// Unit tests for loop_deform (computeU / applyLoopDeformation)
//
// Build & run:
//   g++ -O2 -std=c++14 -I.. -o test_loop_deform test_loop_deform.cpp \
//       ../loop_deform.cpp && ./test_loop_deform
//
// Or via the Makefile:
//   make test
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include "../loop_deform.hpp"

// Helper: max entry-wise asymmetry of M
static double asymmetry(const Mat3& M)
{
    double a = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a = std::max(a, std::abs(M[i][j] - M[j][i]));
    return a;
}

// Helper: Frobenius norm of (M - I)
static double distFromIdentity(const Mat3& M)
{
    double s = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double diff = M[i][j] - (i == j ? 1.0 : 0.0);
            s += diff * diff;
        }
    return std::sqrt(s);
}

// Helper: Frobenius norm of (A - B)
static double frobDiff(const Mat3& A, const Mat3& B)
{
    double s = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double d = A[i][j] - B[i][j];
            s += d * d;
        }
    return std::sqrt(s);
}

// Helper: 3x3 determinant
static double det3(const Mat3& M)
{
    return  M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
}

// =============================================================================
// Test 1: single loop with b = a·[1,0,0], n = [0,1,0], A = L²
//
//   β^p = (a·L²/V) · outer([1,0,0],[0,1,0]) = [[0, a/L, 0],[0,0,0],[0,0,0]]
//   F   = I + β  →  simple shear in the xz-plane
//
//   Expected: U is symmetric (rotation removed), and the (0,1) shear
//   component equals   γ/2 = (a/(2L)) + O(γ³).
// =============================================================================
static void test1_known_result()
{
    const double a = 3.3, L = 100.0;
    const double V = L*L*L;
    const double gamma = a / L;   // β[0][1] = a*L²/L³ = a/L

    LoopDesc lp;
    lp.b[0] = a; lp.b[1] = 0.0; lp.b[2] = 0.0;
    lp.n[0] = 0.0; lp.n[1] = 1.0; lp.n[2] = 0.0;
    lp.A = L * L;

    Mat3 U = computeU({lp}, V);

    // U must be symmetric
    double asym = asymmetry(U);
    assert(asym < 1e-14 && "Test 1 FAILED: U is not symmetric");

    // Off-diagonal shear component: U[0][1] = γ/2 to first order
    double expected_shear = gamma / 2.0;
    double actual_shear   = U[0][1];
    // For γ = 0.033, the second-order correction is O(γ³) < 1e-5
    assert(std::abs(actual_shear - expected_shear) < 1e-5 &&
           "Test 1 FAILED: shear component deviates from γ/2");

    // z-diagonal must be exactly 1 (the z-eigenvector of C is (0,0,1))
    assert(std::abs(U[2][2] - 1.0) < 1e-14 && "Test 1 FAILED: z-diagonal != 1");
    // Off-diagonal elements not involving z must be zero (no coupling to z)
    assert(std::abs(U[0][2]) < 1e-14 && "Test 1 FAILED: U[0][2] != 0");
    assert(std::abs(U[1][2]) < 1e-14 && "Test 1 FAILED: U[1][2] != 0");

    std::cout << "Test 1 PASSED  [known result]"
              << "  U[0][1] = " << actual_shear
              << "  (expected γ/2 = " << expected_shear << ")"
              << "  asym = " << asym << "\n";
}

// =============================================================================
// Test 2: volume conservation for a pure shear loop (b · n = 0)
//
//   b = a/√2 · [1,1,0],  n = [1,-1,0]/√2   →  b·n = 0
//   det(U) must equal 1 to machine precision.
// =============================================================================
static void test2_volume_conservation()
{
    const double a0 = 3.3024, L = 50.0 * a0;
    const double V  = L * L * L;

    // b = (a0/√2) · [1,1,0]  →  components (a0/√2 / √2, a0/√2 / √2, 0) = (a0/2, a0/2, 0)
    // n = [1,-1,0]/√2
    const double bm = a0 / std::sqrt(2.0);
    LoopDesc lp;
    lp.b[0] = bm / std::sqrt(2.0);   // a0/2
    lp.b[1] = bm / std::sqrt(2.0);   // a0/2
    lp.b[2] = 0.0;
    lp.n[0] =  1.0 / std::sqrt(2.0);
    lp.n[1] = -1.0 / std::sqrt(2.0);
    lp.n[2] = 0.0;
    lp.A    = L * L * 0.3;

    // Verify b · n = 0
    double bn = lp.b[0]*lp.n[0] + lp.b[1]*lp.n[1] + lp.b[2]*lp.n[2];
    assert(std::abs(bn) < 1e-14 && "Test 2 setup: b·n != 0");

    Mat3 U = computeU({lp}, V);
    double detU = det3(U);

    assert(std::abs(detU - 1.0) < 1e-10 &&
           "Test 2 FAILED: det(U) != 1 for pure shear loop");

    std::cout << "Test 2 PASSED  [volume conservation]"
              << "  det(U) = " << std::setprecision(15) << detU << "\n"
              << std::setprecision(6);
}

// =============================================================================
// Test 3: superposition — two identical loops == one loop with double area
//
//   β is accumulated linearly, so F and U must be exactly identical when
//   two loops with area A are provided vs. one loop with area 2A.
// =============================================================================
static void test3_superposition()
{
    const double a0 = 3.3, L = 100.0;
    const double V  = L * L * L;

    LoopDesc base;
    base.b[0] = a0 * 0.5; base.b[1] = a0 * 0.3; base.b[2] = a0 * 0.1;
    base.n[0] = 0.0;       base.n[1] = 0.0;       base.n[2] = 1.0;
    base.A    = L * L * 0.1;

    // Two loops, each with area A
    LoopDesc lp1 = base, lp2 = base;
    Mat3 U_two = computeU({lp1, lp2}, V);

    // One loop with area 2A
    base.A *= 2.0;
    Mat3 U_one = computeU({base}, V);

    double diff = frobDiff(U_two, U_one);
    assert(diff < 1e-14 && "Test 3 FAILED: two loops != one doubled loop");

    std::cout << "Test 3 PASSED  [superposition]"
              << "  ||U_two - U_one||_F = " << diff << "\n";
}

// =============================================================================
// Test 4: empty loop list → U = I, H and positions unchanged
// =============================================================================
static void test4_zero_loops()
{
    std::vector<LoopDesc> loops;   // empty

    // Box: arbitrary non-cubic
    Mat3 H = {};
    H[0][0] = 50.0; H[1][1] = 60.0; H[2][2] = 70.0;
    const Mat3 H_orig = H;

    // Positions: a few atoms
    std::vector<double> ax = {1.0,  5.0, -3.0};
    std::vector<double> ay = {2.0, -4.0,  7.0};
    std::vector<double> az = {3.0,  6.0,  1.0};
    const auto ax0 = ax, ay0 = ay, az0 = az;

    applyLoopDeformation(loops, H, ax, ay, az, 3);

    // H must be unchanged
    assert(frobDiff(H, H_orig) < 1e-14 && "Test 4 FAILED: H changed");

    // Positions must be unchanged
    for (int k = 0; k < 3; k++) {
        assert(std::abs(ax[k] - ax0[k]) < 1e-14);
        assert(std::abs(ay[k] - ay0[k]) < 1e-14);
        assert(std::abs(az[k] - az0[k]) < 1e-14);
    }

    std::cout << "Test 4 PASSED  [zero loops → identity]\n";
}

// =============================================================================
// Test 5: after polar decomposition, U must be symmetric  |U - U^T| < 1e-12
//         (tested with a general non-shear loop to stress the algebra)
// =============================================================================
static void test5_rotation_removal()
{
    const double a0 = 3.3, L = 50.0;
    const double V  = L * L * L;

    // General loop (b · n ≠ 0 in general — tests the full polar decomposition)
    LoopDesc lp;
    lp.b[0] = a0 * 0.5; lp.b[1] = a0 * 0.3; lp.b[2] = a0 * 0.2;
    lp.n[0] = 0.2;       lp.n[1] = 0.8;       lp.n[2] = 0.3;
    lp.A    = L * L * 0.4;

    Mat3 U = computeU({lp}, V);

    double asym = asymmetry(U);
    assert(asym < 1e-12 && "Test 5 FAILED: U is not symmetric");

    std::cout << "Test 5 PASSED  [rotation removed, U symmetric]"
              << "  max|U - U^T| = " << asym << "\n";
}

// =============================================================================
int main()
{
    std::cout << std::setprecision(6);
    std::cout << "Running loop_deform unit tests...\n\n";

    test1_known_result();
    test2_volume_conservation();
    test3_superposition();
    test4_zero_loops();
    test5_rotation_removal();

    std::cout << "\nAll tests passed.\n";
    return 0;
}
