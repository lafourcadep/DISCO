# CPU or GPU

Both backends evaluate the same pair term (`Code/kroupa.hpp`), in double precision.

**GPU (`exec_cuda`).** One thread per atom. The loopon positions of the current loop are loaded
cooperatively into shared memory. A tiled kernel handles loops too large for shared memory.

**CPU (`exec_cpu`, or `exec_cuda --cpu`).** The atoms are shared among the OpenMP threads. For each
atom and periodic image, the loop over the loopons is vectorised (SIMD). The summation order is
the GPU's, up to the SIMD partial sums.

## Same sample

Five test cases were run on both backends: a single loop, a parallelogram, a hexagonal prismatic dipole,
a tube quadrupole, and one 6705-loopon loop on the GPU tiled path. The CPU `.xyz` files are
byte-identical to the GPU ones, except for one atom in one case that differs in the last printed
digit (1e-8 Å).

## Speed

Speed is measured in atom–loopon pairs per second, where each pair includes its 27 periodic images.
The run time is roughly atoms × total loopons / rate.

| Hardware | Backend | Rate |
|---|---|---|
| NVIDIA H100 | GPU | 6.1 G pairs/s |
| Intel Core Ultra 5 235H, 14 threads (AVX2) | CPU | 83 M pairs/s (7.7 M on 1 thread) |
| NVIDIA RTX PRO 1000 Blackwell laptop GPU | GPU | 51 M pairs/s |

Data-centre GPUs (A100, H100) are fast in double precision. Workstation and laptop GPUs are not,
and a multi-core CPU can beat them (last row). The CPU backend scales almost linearly with the
number of cores.

!!! tip "Timing runs"
    Use `Output_format: none` to time the displacement computation without writing a file.
