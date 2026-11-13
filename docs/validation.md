# Validation

Every check that DISCO has passed, in one place. They fall into four groups: analytic unit tests
(`make test`), cross-checks between the two field methods, cross-checks between backends and
precisions, and atomistic analysis (OVITO DXA, CNA, PTM) of generated samples. The ground truth of the
generated samples is the JSON sidecar written by the [generators](generators.md).

## Analytic unit tests

`make test` (in `Code/`) builds and runs seven CPU-only tests.

**`tests/test_barnett.cpp`**: the line-integral field (`Method: line`).

| # | Test | Check |
|---|---|---|
| 1 | Against a 7-point Gauss quadrature of the Kroupa kernel over the cut | < 10⁻⁶ Å (measured ≤ 4 × 10⁻¹⁰ Å), glide, prismatic and mixed loops, every shape |
| 2 | Jump across the cut | −**b** inside the contour, 0 outside |
| 3 | Invariance | fan apex moved, starting vertex rotated: < 10⁻¹¹ Å |
| 4 | Loopon sum | converges to the line integral as `Loopon_size` → 0 |
| 5 | Far field | a small loop matches one loopon of the same area |
| 6 | Polyline | counter-clockwise about **n**, area of the contour |

The sign convention of the line integral was fixed by test 1: of the 16 sign combinations of its three
terms and the contour orientation, only one agrees with the quadrature (to 2 × 10⁻¹⁰ Å).

**`tests/test_prismatic.cpp`**: loops with **b** · **n** ≠ 0. The Kroupa jump across the cut is
−**b** · n̂, so **m** · **n** > 0 is a vacancy loop and **m** · **n** < 0 an interstitial one; the slab holds
|**b** · n̂|/d atomic planes and ρ A |**b** · n̂| atoms. Covered: BCC ½⟨111⟩ on (111) and on {110},
a⟨100⟩ on (100); FCC ½⟨110⟩ on (110) and on {111}, a⟨100⟩ on (100).

**`tests/test_snap.cpp`**: by brute force over the generated BCC and FCC sites, the loop centre is
moved exactly halfway between two atomic planes, for {110}, {112}, {123}, {100}, {111} and {113}
normals, with even and odd `Nx`, `Ny`, `Nz`.

**`tests/test_orient.cpp`**: oriented crystals. The atom count is atoms per cell × V / a0³ and every
atom has the perfect first-neighbour shell (12 FCC, 8 BCC) across the periodic faces (the box is
commensurate); the oriented mid-plane snap puts every site at a half-integer number of spacings;
the crystal ↔ box rotations are inverse.

**`tests/test_polycrystal.cpp`**: one grain with a commensurate orientation reproduces the oriented
crystal site for site; with random grains every atom lies in its grain's Voronoi cell, no pair is
closer than the GB cutoff, the build is deterministic; orientations are proper rotations and Euler
angles round-trip; the portable random generator reproduces mt19937_64; loop clipping keeps an inside
polygon, turns an oversized one into the cell's cross-section (area checked on a grid), and empties a
polygon away from the cell.

**`tests/test_twin.cpp`**: nearest lattice site; FCC and BCC lamellae of 3, 6, 9 layers keep the atom
count and leave the matrix unchanged; FCC atoms keep 12 neighbours at the perfect distance; every lamella
atom is the mirror image of a matrix atom across both TBs; thickness snapped to a multiple of 3 layers;
twin frame = parent turned by 180° about **n**.

**`tests/test_loop_deform.cpp`**: the (unused) box deformation: simple shear gives U symmetric with
U₀₁ = γ/2; det U = 1 for **b** · **n** = 0; two loops of area A ≡ one of area 2A; no loop gives U = I;
max |U − Uᵀ| < 10⁻¹² for a general loop.

## Line integral against loopons

On a glide quadrupole, the loopon field converges to the line integral as the grid is refined:
rms difference 9 × 10⁻³ Å at `Loopon_size` = 0.4 a0 and 2.6 × 10⁻⁴ Å at 0.05 a0, largest at the cores.

## Backends and precision

- **CPU and GPU, double precision.** Five cases (single loop, parallelogram, hexagonal prismatic
  dipole, tube quadrupole, one 6705-loopon loop on the GPU tiled path): the CPU `.xyz` files are
  byte-identical to the GPU ones, except one atom in one case that differs in the last printed digit
  (10⁻⁸ Å).
- **Single precision (`--fp32`) against double.** `Method: line`: differences below 3 × 10⁻⁵ Å
  (rms ~10⁻⁶ Å) on samples of up to 315 000 atoms. `Method: loopons`: at most 5 × 10⁻⁶ Å
  (rms ~2 × 10⁻⁷ Å) on four samples of up to 128 649 atoms, including interstitial loops with atoms
  d/2 from the cut; DXA finds the same loops. Polycrystal with twin-aware loops: GPU fp32 = CPU within
  4.3 × 10⁻⁶ Å.
- **GPU chunking** leaves each atom's sum unchanged: identical output with and without it.

## Generated samples (OVITO)

All analyses on unrelaxed samples written by DISCO, unless stated.

**Glide loops and quadrupoles (BCC).** DXA finds each loop of every shape with its ½⟨111⟩ Burgers
vector; tube quadrupoles show the a⟨100⟩ junction legs (or ⟨110⟩ products) predicted by the sidecar.

**Prismatic dipoles from shear loops (BCC).** Three hexagonal dipoles (side 25.6 Å, random **b** and
separation): DXA finds exactly 6 closed ½⟨111⟩ loops, two per dipole with ±**b**, at the predicted
positions. Same in a 61 × 51 × 41 box with `--margin 6` (odd cell counts, different per axis, no loop
crossing the boundary).

**Prismatic loops with matter (BCC).** Interstitial and vacancy hexagons on (111) (R = 22 Å): adaptive
CNA finds no extra non-BCC atoms at the corners and DXA finds one closed ½⟨111⟩ loop, with or without
`Corner_radius`. With the generator's default ⟨112⟩ sides (b = ½[111], R = 22 Å, 217 atoms added): DXA
finds one closed ½⟨111⟩ loop of length 131.3 Å (contour perimeter 132 Å).

**FCC loops.** 22 loops from `make_fcc_loops.py` (single rhombus, 8 mixed rhombi, 3 vacancy Frank
triangles, 10 mixed Frank and prismatic loops): one DXA loop per generated loop with its Burgers vector
(½⟨110⟩ or ⅓⟨111⟩) and no other segment, HCP fault layers inside the Frank loops only, and atom counts
within the contour scatter of ρ A |**b** · n̂|. `Examples/fcc_faulted_loops.yaml`: DXA finds ⅓⟨111⟩ Frank and
⅙⟨112⟩ Shockley loops with HCP layers inside.

**Oriented crystals.** FCC (X = [1−10], Y = [11−2], Z = [111]) and BCC (X = [111], Y = [−110],
Z = [−1−12]): PTM finds the requested orientation on every sampled atom of a perfect box, there are no
defects on the periodic faces, and DXA finds each loop with its Burgers vector rotated into the box frame.

**Polycrystals.** Cu (8 random grains, 120 Å box) and Ta (4 listed grains):

- the `grain` column is the nearest seed for every atom;
- atoms more than 5 Å from any GB are 100 % FCC (99.985 % BCC, the rest at triple lines), so grains
  crossing the box faces have no seam;
- the PTM orientation of the interior atoms of each grain matches its R (400/400 and 200/200 sampled
  atoms); listed orientations and Euler angles are reproduced;
- the minimum interatomic distance is `GB_min_distance` × NN; OVITO's grain segmentation finds the
  expected number of grains.

**Loops in grains.**

- `polycrystal_cu_loops.yaml`: atoms are added or removed only in the grains of the Frank and prismatic
  loops; extra HCP atoms appear only in the grains of the Shockley (GB-to-GB fault) and Frank loops; DXA
  finds the glide, Frank and prismatic loops with Burgers vector R_k **b**. The GB-to-GB Shockley fault
  has no DXA line: all its partials lie in GBs, which DXA does not analyse.
- Generated Cu sample (12 grains, 41 loops, 4.7 × 10¹⁶ m⁻² after clipping, 658 k atoms): DXA finds all 41
  loops in their grains with Burgers vector R_k **b** (within 8°: the loops' own elastic lattice rotation
  reaches a few degrees near the cores), plus short segments where clipped Frank loops meet the
  unrelaxed GBs.
- Generated Ta sample (6 grains, 8 loops): atoms added or removed only in the loops' grains; DXA finds
  the loops away from GBs. On *unrelaxed* BCC GBs DXA reports short spurious segments even without
  loops and can miss a loop within ~15 Å of such a GB; there the displacement jump across the cut still
  shows the loop (−2.4 Å along **b** for |**b**| = 2.86 Å, at the loop edge). Relax before analysis.

**Twins.**

- FCC single crystal (two lamellae, m = 6 and 9): exactly one HCP layer per TB (2688 atoms), all other
  atoms FCC, lamella orientation = R_twin within 0.012°, 60.00° from the matrix, no DXA line. BCC (112),
  m = 9: lamella 60.00° from the matrix, no DXA line (three non-BCC layers per unrelaxed reflection TB).
- Polycrystals (3 hand-placed and 7 random lamellae): every lamella at 60.00° from its grain, atoms
  outside lamellae unmoved, largest shift a0/√6. Where a lamella ends on an unrelaxed GB, DXA finds a few
  short junction segments (about 36 Å in total for 3 lamellae).
- Twin-aware loops: a glide and a Shockley loop inside a lamella are found with **b** = R_twin **b**_twin
  (0.7° and 0.7°), a matrix loop clipped at a TB with **b** = R **b** (0.0°); in a polycrystal lamella
  0.3° from R_twin **b** (70.6° from the grain-frame **b**).

## Not yet done

- MD stability at finite temperature (relaxed box, 600 K) of prismatic and Frank loops, counted with
  CNA / DXA.
- `exec_cuda` on an NVIDIA GPU since the line-integral method: fp64, `--fp32` and `exec_cpu` compared
  for both methods.
