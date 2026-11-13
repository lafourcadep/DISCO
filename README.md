# DISCO — DISlocation COnfiguration generator

[![Documentation](https://github.com/lafourcadep/DISCO/actions/workflows/docs.yml/badge.svg)](https://lafourcadep.github.io/DISCO/)

DISCO builds atomistic BCC or FCC samples containing dislocation loops, for molecular dynamics.
The displacement of every atom is the isotropic elastic field of the loops, with periodic images:
by default an exact closed-form line integral along each contour, or a sum over infinitesimal loops
("loopons", Kroupa 1966) tiling each loop. It runs on a GPU (CUDA, or Apple GPUs through mvcc) or
on the cores of a CPU (OpenMP). Loops can be circles, rectangles, parallelograms or convex polygons
(hexagons). Prismatic loops (interstitial or vacancy) get the matching atoms added or removed.
Oriented crystals, Voronoi polycrystals (loops clipped at the grain boundaries) and coherent twin
lamellae are supported. Python generators write inputs, each with a JSON ground truth: BCC single
loops on any slip system, junction and annihilation quadrupoles, prismatic dipoles and prismatic
loops; FCC prismatic and Frank loops; populated FCC or BCC polycrystals.

**Documentation: <https://lafourcadep.github.io/DISCO/>**

## Quick start

Requirements: a C++14 compiler with OpenMP, [yaml-cpp](https://github.com/jbeder/yaml-cpp);
the CUDA toolkit for the GPU build; Python 3 + numpy for the generators.

```bash
cd Code
make exec_cpu                                   # CPU build
make GPU_ARCH=sm_80                             # GPU build (exec_cuda), set your architecture
make test                                       # unit tests
./exec_cpu ../Examples/quadrupole_circle.yaml   # writes Geo_40x40x40_with_loop.lmp
```

## Layout

| Folder | Contents |
|---|---|
| `Code/` | C++/CUDA sources, Makefile, unit tests, `analysis_dxa.py` (OVITO DXA) |
| `Examples/` | Input templates, `generators/` (input generators with ground-truth sidecars), a LAMMPS input template |
| `Parameters/BCC.txt` | BCC slip systems used by the generators (same numbering as the DANCE analysis code) |
| `docs/` | Documentation sources ([Zensical](https://zensical.org)); built and published by GitHub Actions |

Build the documentation locally with `pip install zensical && zensical serve`.

## License

MIT, see [LICENSE](LICENSE).
