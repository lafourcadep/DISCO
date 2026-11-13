# Installation

## Requirements

| | Needed for |
|---|---|
| C++14 compiler with OpenMP (g++) | both executables |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp) | input parser, both executables |
| CUDA toolkit (`nvcc`) | `exec_cuda` only |
| [mvcc](https://github.com/doximity/mvcc) toolkit, Homebrew `libomp` | `exec_mvcc` only (macOS, Apple silicon) |
| Python 3 + numpy | the input generators in `Examples/generators/` |

## Build

All build commands run in `Code/`:

```bash
make            # ./exec_cuda: GPU build (needs nvcc; set GPU_ARCH, see below)
make exec_cpu   # ./exec_cpu: CPU-only build (g++ with OpenMP, no CUDA toolkit)
make all        # both
make exec_mvcc  # ./exec_mvcc: Apple-silicon GPU build through mvcc (see below)
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

### Apple silicon GPU (mvcc)

[mvcc](https://github.com/doximity/mvcc) compiles CUDA for Apple GPUs (Metal), so `exec_mvcc` is
built from the same `backend_cuda.cu`. Build mvcc's toolkit (`tools/install_toolkit.sh` in its
repository), then:

```bash
make exec_mvcc MVCC_NVCC=/path/to/mvcc/toolkit/bin/nvcc   # or put toolkit/bin first in PATH
./exec_mvcc --fp32 ../Examples/quadrupole_circle.yaml
```

Apple GPUs have no fp64 hardware: without `--fp32` the kernels run in emulated double precision,
far slower than `exec_cpu`. Long kernels are split automatically into launches of a few seconds
(see [CPU or GPU](backends.md#command-line-options)). Tested with mvcc 1.1.1 on an M4 Pro.

## Unit tests

```bash
make test
```

Seven tests: the line-integral field against a quadrature of the Kroupa kernel (`test_barnett`), the
matter slab of prismatic loops (`test_prismatic`), the mid-plane snap (`test_snap`), oriented crystals
(`test_orient`), polycrystals and loop clipping (`test_polycrystal`), twin lamellae (`test_twin`) and
the box-deformation routine (`test_loop_deform`). What each one checks:
[Validation](validation.md#analytic-unit-tests).
