# DISCO

**DISlocation COnfiguration generator.** DISCO builds atomistic BCC samples that contain
dislocation loops, ready for molecular dynamics.

The displacement of every atom is the isotropic elastic field of all loops, with periodic images.
By default it is computed exactly, as a closed-form line integral along each loop contour. It can
also be computed as the sum over small square (or rhombic) elements tiling each loop, the *loopons*
(Kroupa, 1966). The computation runs on a GPU (CUDA, or Apple GPUs through mvcc) or on the cores of a
CPU (OpenMP).

<figure markdown="span">
  ![A rounded-rectangle loop tiled by loopons](assets/kroupa_loopons.png){ width="380" }
  <figcaption>A dislocation loop for the slip system [−111](110), tiled by loopons (black grid).
  m is the Burgers vector direction, n the loop-plane normal.</figcaption>
</figure>

## What it does

- BCC lattice of `Nx × Ny × Nz` cubic cells, any lattice parameter.
- Any number of loops: circles, rectangles (optionally rounded, optionally rotated in their plane),
  parallelograms and convex polygons such as hexagons (optionally rounded), on any plane `n` with
  Burgers vector `b = a0/2 m`.
- Prismatic and mixed loops (b · n ≠ 0), interstitial or vacancy: the atoms of the slab inside the
  loop are added or removed, so the loops survive thermal MD.
- Loop centres are moved onto the mid-plane between two atomic planes of the loop plane family.
- Output as extended XYZ (OVITO) or LAMMPS data file.
- Python generators for common configurations: single loops on all 48 BCC slip systems,
  junction and annihilation quadrupoles (one, or many different ones in one sample),
  prismatic dipoles (one, or N at random positions), prismatic loops with added or removed atoms
  (circular or hexagonal), and large tiled samples with random positions. Each generated input comes with a JSON file holding the expected
  dislocation content.

| | | | |
|:-:|:-:|:-:|:-:|
| ![](assets/g1_circle.png){ width="110" } | ![](assets/g2_rectangle.png){ width="110" } | ![](assets/g3_parallelogram.png){ width="110" } | ![](assets/g4_rect112.png){ width="110" } |
| {110} circle | {110} rounded rectangle | ⟨111⟩ parallelogram | {112} rectangle |
| ![](assets/g5_tube_junction.png){ width="110" } | ![](assets/g6_tube_annihilation.png){ width="110" } | ![](assets/g7_prismatic.png){ width="110" } | ![](assets/g8_hexagonal.png){ width="110" } |
| tube, junction | tube, annihilation | prismatic dipole | hexagonal interstitial loop |

*Dislocation lines extracted by OVITO's DXA from samples written by DISCO.*

## Workflow

```mermaid
flowchart LR
  G["Python generators<br/>(YAML input)"] --> D["DISCO<br/>exec_cuda / exec_cpu"]
  D -->|".lmp / .xyz"| L["LAMMPS<br/>relaxation, loading"]
  L -->|"dumps"| O["OVITO DXA<br/>analysis_dxa.py"]
  O -->|".ca"| A["Network analysis<br/>(e.g. DANCE)"]
```

## Quick start

```bash
git clone https://github.com/lafourcadep/DISCO.git
cd DISCO/Code
make exec_cpu                                   # CPU build (g++ + OpenMP + yaml-cpp)
./exec_cpu ../Examples/quadrupole_circle.yaml   # writes Geo_40x40x40_with_loop.lmp
```

With a CUDA toolkit, `make GPU_ARCH=sm_80` builds `exec_cuda` (see [Installation](installation.md)).

<figure markdown="span">
  ![Dislocation network](assets/network_32M.png){ width="460" }
  <figcaption>Early state of the dislocation network in a 32-million-atom BCC Ta sample
  seeded with randomly placed dislocation quadrupoles.</figcaption>
</figure>

## License

MIT. See [`LICENSE`](https://github.com/lafourcadep/DISCO/blob/main/LICENSE).
