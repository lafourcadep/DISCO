# CPU or GPU

Both backends evaluate the same pair term (`Code/kroupa.hpp`), in double precision.
They only receive loopon positions and loop parameters (b, m, n, dS), so they do not depend on
the loop shape. The atoms added or removed for prismatic loops are handled on the host, before
and after the kernel, identically for both backends
(see [Theory](theory.md#prismatic-loops-adding-or-removing-matter)).

**GPU (`exec_cuda`).** One thread per atom. The loopon positions of the current loop are loaded
cooperatively into shared memory. A tiled kernel handles loops too large for shared memory.

**Apple-silicon GPU (`exec_mvcc`).** The same CUDA source, compiled with
[mvcc](https://github.com/doximity/mvcc), which translates CUDA to Metal (see
[Installation](installation.md#apple-silicon-gpu-mvcc)). Apple GPUs have no fp64 hardware, so run it
with `--fp32` (below).

**CPU (`exec_cpu`, or `exec_cuda --cpu`).** The atoms are shared among the OpenMP threads. For each
atom and periodic image, the loop over the loopons is vectorised (SIMD). The summation order is
the GPU's, up to the SIMD partial sums.

## Command-line options

```text
exec_cuda [--cpu] [--fp32] [--gpu-chunk <pairs>] <input.yaml>
```

| Option | Effect |
|---|---|
| `--cpu` | OpenMP backend even if a GPU is present |
| `--fp32` | GPU only: single-precision pair term (see below) |
| `--gpu-chunk <pairs>` | GPU only: at most this many atom–loopon pairs per kernel launch, the host waiting after each launch. `0` turns chunking off. Default: on for Apple GPUs (10⁷ pairs in fp64, 2 × 10⁸ in fp32), off otherwise |

**Chunking.** A launch then covers only a range of atoms. Each atom's sum is unchanged, so the
output is identical with or without chunking. It exists because mvcc's runtime aborts any GPU wait
longer than 8 s (reported as `an illegal memory access was encountered`). On NVIDIA GPUs it is
not needed.

**Single precision (`--fp32`).** The pair term is evaluated in float, with three safeguards:

- loopons are uploaded as float offsets from the loopon centroid c;
- each atom is reduced to a − c in double before rounding to float, so coordinates near the loop,
  where the field varies fastest, keep their precision;
- the sum over the loopons of each periodic image is in float, and the 27 image sums in double.

Loopons in float also halve the shared memory, so loops of up to 2730 loopons take the fast path.
Compared with the double-precision output, displacements differ by at most 5 × 10⁻⁶ Å
(rms ~2 × 10⁻⁷ Å). This was measured on four samples of up to 128 649 atoms, including
interstitial loops with atoms d/2 from the cut, and DXA finds the same loops.

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
| Apple M4 Pro, 12 threads | CPU | 182 M pairs/s |
| Apple M4 Pro GPU (16 cores), mvcc 1.1.1, `--fp32` | GPU | 540–810 M pairs/s |
| Apple M4 Pro GPU (16 cores), mvcc 1.1.1, double | GPU | ~4 M pairs/s (fp64 emulated) |

Data-centre GPUs (A100, H100) are fast in double precision. Workstation and laptop GPUs are not,
and a multi-core CPU can beat them (RTX PRO 1000 row). On an Apple GPU only `--fp32` is worthwhile:
it is 3–4.5× faster than the 12 CPU cores, while double precision is ~40× slower. The CPU backend scales almost linearly with the
number of cores.

!!! tip "Timing runs"
    Use `Output_format: none` to time the displacement computation without writing a file.
