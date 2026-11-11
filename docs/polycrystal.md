# Polycrystals

A `polycrystal:` block replaces `Nx`, `Ny`, `Nz` (and `Orientation`) and builds a periodic
Voronoi polycrystal: each grain has a seed and an orientation, and each atom belongs to the grain
of its nearest seed (periodic minimum image).


```yaml
general:
  Lattice: FCC
  Element: Cu
  a0: 3.615
  Output_format: both

polycrystal:
  Box: [200, 200, 200]          # Å; any size (the GBs absorb the lattice mismatch)
  Grains: 12                    # random seeds and orientations ...
  Seed: 7                       # ... from this random seed (default 1)
  GB_min_distance: 0.8          # optional: pairs closer than 0.8 × NN distance lose one atom
```

or an explicit list (each key optional: a missing centre or orientation is random):

```yaml
polycrystal:
  Box: [100, 90, 110]
  Seed: 2
  Grain_list:
    - Center: [-25, -20, -30]                                      # Å, box-centred
      Orientation: {X: [1, 1, 1], Y: [-1, 1, 0], Z: [-1, -1, 2]}   # crystal directions of x, y, z
    - Center: [25, 20, 30]
      Euler: [30, 45, 60]                                          # Bunge φ1 Φ φ2, degrees
    - Center: [25, -25, -20]                                       # random orientation
    - {}                                                           # random centre and orientation
```

## Construction

- **Orientations.** Random orientations are uniform over all rotations (random unit quaternions).
  `Orientation` takes integer directions as in [Oriented crystal](input.md#oriented-crystal);
  `Euler` takes Bunge angles. Each grain's lattice has a site on its seed.
- **Filling.** A grain's periodic Voronoi cell lies within ± L/2 of its seed on every axis, so each
  grain is filled from its own lattice in that region, unwrapped around the seed, and wrapped back
  into the box. A grain that crosses a box face is one continuous crystal: there is no seam at the
  periodic boundary. Ties between two seeds are broken by a fixed rule, so no site is kept twice;
  a single grain with a commensurate orientation reproduces the oriented crystal exactly.
- **Grain boundaries.** Sites of neighbouring grains can come too close. Every pair closer than
  `GB_min_distance` × NN distance (FCC a0/√2, BCC a0 √3/2) loses one atom (the one of the higher
  grain). A perfect crystal has no such pair, so only GB atoms are affected. The log gives the number
  removed; it is a few % for grains of ~10 nm and grows as grains get smaller.
- **Random numbers** are platform independent: the same `Seed` gives the same sample with any
  compiler.

!!! warning "Relax before use"
    The GBs are built geometrically, not relaxed. Minimise the energy (and anneal if needed) before
    any MD. With few grains a grain can meet its own periodic image, which forms a GB like any other.

## Loops in grains

Each item of `dislocations` takes `Grain: k` (1-based, as in the `grain` column):

```yaml
dislocations:
  - Grain: 4                 # crystal frame of grain 4
    Shape: circular
    Radius: 80.0             # larger than the grain: see below
    Center: [43.985, -27.887, 14.467]   # box frame, origin at the box centre
    m: [-1/3, -1/3, 2/3]     # Shockley loop, stable sense for n = [111]
    n: [1, 1, 1]
    Clip: auto               # optional: auto (default) | none
```

- `m`, `n`, `Side`, `Dir_1`, `Dir_2`, `Vertex_dir`, `Vertices` are in the **crystal frame of grain k**;
  `Center` is in the box frame. The centre is snapped to a mid-plane of grain k's lattice (from the
  seed, so it also works next to a box face).
- **Clipping.** The loop's cut surface must stay in its grain, where b is a lattice vector. If the
  contour leaves the grain, it is clipped by the grain's Voronoi faces: the loop keeps the part of its
  contour inside the grain, and the clipped edges lie in the GBs (dislocation segments absorbed by
  the boundary). A loop larger than its grain therefore spans the whole grain cross-section: a
  Shockley loop becomes a stacking fault from GB to GB, a perfect or Frank loop cuts the grain
  through. A clipped loop is rebuilt as a polygon (circles use their chords, `Segment_tolerance`).
- **Skipped loops.** A loop whose plane misses its grain (or keeps a sliver smaller than 4 |b|²) is
  skipped, with a warning. Nothing is lost: a loop's whole displacement field, rotations included,
  comes from its dislocation lines, and an empty contour has none; keeping it unclipped would put its
  cut, and a fault, in another grain.
- **Field.** Each loop's field is applied to every atom of every grain: with isotropic elasticity the
  field does not depend on the crystal orientation, so the polycrystal is elastically homogeneous.
- **Atoms added or removed** (b · n ≠ 0, `Matter: auto`) are taken in the loop's grain only.
- `Clip: none` keeps the contour as given (warning if it leaves the grain): its cut crosses the GB.

Example: `Examples/polycrystal_cu_loops.yaml` (one loop of each case above).
For populations, see `make_polycrystal_loops.py` in [Input generators](generators.md).

## Twin lamellae

A top-level `twins:` list adds coherent twin lamellae, in grains or in an oriented single crystal:

```yaml
twins:
  - Grain: 3                 # polycrystal: the grain (1-based); single crystal: omit
    n: [1, 1, 1]             # twin plane, crystal frame of the grain: {111} FCC, {112} BCC
    Center: [-12.3, -23.0, 39.9]   # box frame: a point on the first twin boundary (TB)
    Thickness: 10.0          # Angstrom, snapped to a multiple of 3 planes (here 12.52 A)
```

- **Construction.** A coherent twin is the mirror image of the parent across a twin plane through
  atoms. The first TB is moved onto the nearest atomic plane; the atoms of the next m − 1 layers are
  moved, each layer by one in-plane vector (at most a0/√6 in FCC), onto the mirror lattice. No atom
  is added or removed and there is no elastic field. The second TB (layer m) is also a coherent
  mirror only if the two mirrors compose to a lattice translation: m must be a multiple of 3 for
  FCC {111} and BCC {112}, and the thickness is snapped accordingly (minimum 3 planes).
- **Extent.** A lamella spans its whole grain and ends at the GBs. In a single crystal it spans the
  periodic box, so the twin plane must be normal to a box axis (use `Orientation`, e.g. Z = [111]).
  Lamellae of one grain must be parallel and separated by at least one matrix layer.
- **Twin frame.** The twin's crystal frame is the parent's turned by 180° about n (Σ3: 60° about
  ⟨111⟩). A loop whose centre lies inside a lamella is read in the **twin frame** (`m`, `n`, `Side`, …
  are twin-crystal directions) and snapped to the twin's planes. Every loop of a twinned grain (or
  crystal) is clipped at the TBs like at the GBs, so its cut stays in its own crystal, and its
  matter slab only takes atoms of its own region.
- **BCC.** The pure-mirror {112} TB (reflection twin) puts atoms of layers ±1 at 2d = 0.816 a0 from
  each other; relaxation turns it into the isosceles TB (shift a0/12 ⟨111⟩). Relax before MD.
- Output: an xyz column `twin` (0 = matrix, t = lamella t), `<name>_twins.txt` (plane, normal,
  TB point, m, thickness, atoms moved, R_twin), and a `region` column in `Poly_..._loops.txt`.

Examples: `Examples/fcc_twin_lamella.yaml` (single crystal, with loops in and around the lamellae),
`Examples/polycrystal_cu_twins.yaml`. For random lamellae in a polycrystal: `make_polycrystal_loops.py
--twins-per-grain` ([Input generators](generators.md)).

Checked with OVITO:

- FCC single crystal (two lamellae, m = 6 and 9): exactly one HCP layer per TB (2688 atoms), all other
  atoms FCC, lamella orientation = R_twin within 0.012°, 60.00° from the matrix, no DXA line. BCC (112),
  m = 9: lamella 60.00° from the matrix, no DXA line (three non-BCC layers per unrelaxed reflection TB).
- Polycrystals (3 hand-placed and 7 random lamellae): every lamella at 60.00° from its grain, atoms
  outside lamellae unmoved, largest shift a0/√6. Where a lamella ends on an unrelaxed GB, DXA finds a
  few short junction segments (about 36 Å in total for 3 lamellae).
- Twin-aware loops: a glide and a Shockley loop inside a lamella are found with b = R_twin b_twin
  (0.7° and 0.7°), a matrix loop clipped at a TB with b = R b (0.0°); in a polycrystal lamella
  0.3° from R_twin b (70.6° from the grain-frame b). GPU fp32 = CPU within 4.3e-6 Å.

## Output

- `Poly_<G>g_<Lx>x<Ly>x<Lz>.xyz` / `.lmp` (box lengths rounded to Å). The extended XYZ file has an
  extra per-atom column `grain` (1-based), which OVITO reads as a particle property. The LAMMPS file
  keeps `atom_style atomic` (one atom type); use the xyz file, or the table below, for grain IDs.
- `Poly_..._grains.txt`: one line per grain: id, seed (box-centred, Å), Bunge Euler angles (degrees),
  the 9 entries of R (rows: crystal directions of the lab x, y, z axes; lab = R · crystal) and the
  number of atoms.
- `Poly_..._loops.txt` (when there are loops): one line per input loop: grain, kind (perfect,
  Frank, Shockley), b in the crystal frame (units of a0) and in the box frame (Å), n, status
  (`kept`, `clipped`, `skipped_outside`, `skipped_sliver`), contour area before and after clipping,
  perimeter, atoms added (+) or removed (−). The log gives the dislocation density (line length of the
  kept loops / volume).
- The log prints the same tables.

## Checks

Checked with OVITO on Cu (8 random grains, 120 Å box) and Ta (4 listed grains):

- the `grain` column is the nearest seed for every atom;
- atoms more than 5 Å from any GB are 100 % FCC (99.985 % BCC, the rest at triple lines), so the
  grains crossing the box faces have no seam;
- the PTM orientation of the interior atoms of each grain matches its R (400/400 and 200/200 sampled
  atoms); listed orientations and Euler angles are reproduced;
- the minimum interatomic distance is `GB_min_distance` × NN;
- OVITO's grain segmentation finds the expected number of grains.

Building 2.2 M Cu atoms in 50 grains takes 0.2 s on 12 CPU threads.

Loops in grains, checked with OVITO:

- `polycrystal_cu_loops.yaml`: atoms are added or removed only in the grains of the Frank and prismatic
  loops; extra HCP (fault) atoms appear only in the grains of the Shockley (end-to-end fault) and
  Frank loops; DXA finds the glide, Frank and prismatic loops with Burgers vector R_k b. The
  end-to-end Shockley fault has no DXA line: all its partials lie in GBs, which DXA does not analyse.
- A generated Cu sample (12 grains, 41 loops, 4.7 × 10¹⁶ m⁻² after clipping, 658 k atoms): DXA finds
  all 41 loops in their grains with Burgers vector R_k b (within 8°: the loops' own elastic lattice
  rotation reaches a few degrees near the cores), plus short segments where clipped Frank loops meet
  the unrelaxed GBs. The field took 55 s on 12 CPU threads.
- A generated Ta sample (6 grains, 8 loops): atoms are added or removed only in the loops' grains, and
  DXA finds the loops away from GBs. DXA on *unrelaxed* BCC GBs reports many short spurious segments
  even without loops, and can miss a loop lying within ~15 Å of such a GB; there the displacement
  jump across the cut plane still shows the loop (−2.4 Å along b for |b| = 2.86 Å, at the loop edge).
  Relax the sample before a dislocation analysis.
