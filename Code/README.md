# Code

Sources of `exec_cuda` (GPU + CPU) and `exec_cpu` (CPU only).

| File | Role |
|---|---|
| `main.cpp` | Driver: input, BCC lattice, loopons, backend choice, output |
| `kroupa.hpp` | Kroupa pair term, shared by both backends |
| `backend.hpp`, `backend_cuda.cu`, `backend_cpu.cpp` | Backend interface, CUDA kernels, OpenMP backend |
| `input.hpp/.cpp` | YAML input parser |
| `loop.hpp/.cpp` | Loop geometry: loopon grid, Burgers vector, mid-plane snap |
| `loop_deform.hpp/.cpp` | Box deformation from the plastic distortion |
| `tests/` | Unit tests (`make test`) |
| `generate_quadrupole.py` | Quadrupole input generator |
| `analysis_dxa.py` | OVITO DXA on LAMMPS dumps (`.ca` / `.vtk` output) |

Build: `make exec_cpu`, `make GPU_ARCH=sm_XX`, `make test`.
Full documentation: <https://lafourcadep.github.io/DISCO/> (sources in `../docs/`).
