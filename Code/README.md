# Code

Sources of `exec_cuda` (GPU + CPU) and `exec_cpu` (CPU only).

| File | Role |
|---|---|
| `main.cpp` | Driver: input, lattice (cubic, oriented, polycrystal), twins, loops, atoms added / removed, backend choice, output |
| `lattice.hpp` | BCC / FCC lattices, oriented crystals, plane spacings, mid-plane snap, Frank / Shockley rules |
| `polycrystal.hpp/.cpp` | Periodic Voronoi polycrystal, GB pruning, loop clipping at the grain faces |
| `twin.hpp/.cpp` | Coherent twin lamellae (mirror shifts, twin frame) |
| `kroupa.hpp` | Kroupa pair term (`Method: loopons`), shared by both backends |
| `barnett.hpp` | Closed-form line-integral field of a polygonal loop (`Method: line`), shared by both backends |
| `backend.hpp`, `backend_cuda.cu`, `backend_cpu.cpp` | Backend interface, CUDA kernels, OpenMP backend |
| `input.hpp/.cpp` | YAML input parser |
| `loop.hpp/.cpp` | Loop geometry: contours (circle, rectangle, parallelogram, convex polygon), polyline, loopon grid, Burgers vector, mid-plane snap |
| `loop_deform.hpp/.cpp` | Box deformation from the plastic distortion |
| `prismatic.hpp/.cpp` | Atoms added / removed for loops with b·n ≠ 0 (interstitial / vacancy) |
| `tests/` | Unit tests (`make test`) |
| `analysis_dxa.py` | OVITO DXA on LAMMPS dumps (`.ca` / `.vtk` output) |

Build: `make exec_cpu`, `make GPU_ARCH=sm_XX`, `make exec_mvcc` (Apple GPU through mvcc), `make test`.
Options: `--cpu`, `--fp32` (GPU, single precision), `--gpu-chunk <pairs>` (GPU launch size).
Full documentation: <https://lafourcadep.github.io/DISCO/> (sources in `../docs/`).
