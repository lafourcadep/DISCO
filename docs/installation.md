# Installation

## Requirements

| | Needed for |
|---|---|
| C++14 compiler with OpenMP (g++) | both executables |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp) | input parser, both executables |
| CUDA toolkit (`nvcc`) | `exec_cuda` only |
| Python 3 + numpy | the input generators in `Examples/generators/` |

## Build

All build commands run in `Code/`:

```bash
make            # ./exec_cuda: GPU build (needs nvcc; set GPU_ARCH, see below)
make exec_cpu   # ./exec_cpu: CPU-only build (g++ with OpenMP, no CUDA toolkit)
make all        # both
make test       # builds and runs the unit tests (CPU only, no GPU needed)
make clean
```

`exec_cpu` is compiled with `-march=native`. Build it on the machine (or the node type)
that will run it.

### Finding yaml-cpp

The Makefile takes the first of these that works:

1. `make YAML_PREFIX=/path/to/yaml-cpp`, which looks in `include/` and in `lib64/`, `lib/`
   or the prefix itself. Alternatively, pass `YAML_INCDIR=…` and `YAML_LIBDIR=…`.
2. An environment module that sets `YAMLCPP_INCDIR` / `YAMLCPP_LIBDIR`
   (e.g. `module load yaml-cpp` on TGCC machines).
3. A spack install of `yaml-cpp@0.8.0` (`spack location -i`).
4. `pkg-config yaml-cpp`.

With options 1–3 the library path is written into the executable (rpath), so no
`LD_LIBRARY_PATH` is needed at run time.

### GPU architecture

Set `GPU_ARCH` to match the GPU, e.g. `make GPU_ARCH=sm_90`:

| GPU | `GPU_ARCH` |
|---|---|
| GTX 1080 / 1070 (Pascal) | `sm_61` |
| V100 (Volta) | `sm_70` |
| RTX 2080 / Quadro RTX (Turing) | `sm_75` |
| A100 (Ampere) | `sm_80` |
| RTX 3090 / 3080 (Ampere) | `sm_86` |
| RTX 4090 / 4080 (Ada Lovelace) | `sm_89` |
| H100 (Hopper) | `sm_90` |
| RTX PRO Blackwell | `sm_120` (Makefile default) |

`nvcc --list-gpu-arch` lists every option your toolkit supports.

## Unit tests

```bash
make test
```

`tests/test_loop_deform.cpp` checks the box-deformation routine (see [Theory](theory.md#box-deformation)):

| # | Test | Check |
|---|---|---|
| 1 | Known result (simple shear) | U symmetric; U[0][1] = γ/2; z block = identity |
| 2 | Volume conservation | det(U) = 1 for b · n = 0 |
| 3 | Superposition | two loops of area A ≡ one loop of area 2A |
| 4 | Empty loop list | U = I; box and positions unchanged |
| 5 | Rotation removal | max\|U − Uᵀ\| < 1e-12 for a general loop |

`tests/test_snap.cpp` checks by brute force, over the generated BCC sites, that the loop centre
is moved exactly halfway between two atomic planes. It covers {110}, {112}, {123} and {100}
normals with even and odd `Nx`, `Ny`, `Nz`.
