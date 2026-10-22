// =============================================================================
// Kroupa dislocation-loop displacement field — CPU (OpenMP) backend
//
// Same work and summation order as the CUDA kernel (backend_cuda.cu): for each
// atom, loop by loop, the 27 periodic images and all loopons, accumulated in
// registers then added to the atom's displacement.
//   • Atoms are distributed over the OpenMP threads (no race: each thread owns
//     its atoms); set the thread count with OMP_NUM_THREADS.
//   • The loopon loop runs over SoA arrays and is vectorised (omp simd); the
//     pair term is kroupaPair (kroupa.hpp), shared with the GPU kernels.
// =============================================================================

#include <iostream>
#include <vector>
#include <omp.h>

#include "backend.hpp"

using namespace std;

int cpuThreads()
{
    return omp_get_max_threads();
}

void kroupaCPU(const AtomArrays& atoms, const vector<LoopData>& loops,
               double dx, double dy, double dz, BackendTiming& t)
{
    const int    nbAtoms = atoms.n;
    const int    nbLoops = (int)loops.size();
    const double k1      = kroupaK1();
    t.loopSec.assign(nbLoops, 0.);

    cout << "Running CPU kernels (" << cpuThreads() << " OpenMP threads)...\n";
    const double tStart = omp_get_wtime();
    for (int l = 0; l < nbLoops; l++) {
        const LoopData& lp  = loops[l];
        const LoopPhys  L   = lp.phys;
        const int       nLp = (int)lp.x.size();
        const double*   lx  = lp.x.data();
        const double*   ly  = lp.y.data();
        const double*   lz  = lp.z.data();
        const double    b_n = L.mx*L.nx + L.my*L.ny + L.mz*L.nz;
        cout << "  Loop " << l+1 << "/" << nbLoops << ":  " << nLp << " loopons\n" << flush;

        const double t0 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 256)
        for (int a = 0; a < nbAtoms; a++) {
            const double ax = atoms.x[a], ay = atoms.y[a], az = atoms.z[a];
            double dux = 0., duy = 0., duz = 0.;

            for (int ii = -1; ii <= 1; ii++) {
            for (int jj = -1; jj <= 1; jj++) {
            for (int kk = -1; kk <= 1; kk++) {
                const double delta_x = dx * ii;
                const double delta_y = dy * jj;
                const double delta_z = dz * kk;

#pragma omp simd reduction(+:dux,duy,duz)
                for (int p = 0; p < nLp; p++)
                    kroupaPair(ax - (lx[p] + delta_x), ay - (ly[p] + delta_y), az - (lz[p] + delta_z),
                               L, k1, b_n, dux, duy, duz);
            }}}

            atoms.ux[a] += dux;
            atoms.uy[a] += duy;
            atoms.uz[a] += duz;
        }
        t.loopSec[l] = omp_get_wtime() - t0;
    }
    t.kernels = omp_get_wtime() - tStart;
    cout << "\nCPU kernel time: " << 1000. * t.kernels << " ms\n\n";
}
