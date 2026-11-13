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

Each grain is filled from its own lattice (with a site on its seed) inside its periodic Voronoi cell;
a grain crossing a box face is one continuous crystal. Pairs of atoms closer than
`GB_min_distance` × NN distance (FCC a0/√2, BCC a0√3/2) lose one atom; only GB atoms are affected, and
the log gives the number removed. Random orientations are uniform over all rotations, and the same
`Seed` gives the same sample on any platform. Details: [Theory](theory.md#polycrystal).

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
  `Center` is in the box frame. The centre is snapped to a mid-plane of grain k's lattice.
- **Clipping** (`Clip: auto`). If the contour leaves the grain, it is clipped by the grain's Voronoi
  faces: the loop keeps the part inside the grain, and the clipped edges lie in the GBs. A loop larger
  than its grain spans the grain's cross-section (a Shockley loop becomes a GB-to-GB stacking fault).
  A clipped loop is rebuilt as a polygon (circles use their chords, `Segment_tolerance`).
- **Skipping.** A loop whose plane misses its grain, or keeps a sliver smaller than 4 |b|², is skipped
  with a warning.
- The field of every loop is applied to every atom (isotropic elasticity: the polycrystal is
  elastically homogeneous). Atoms added or removed (b · n ≠ 0) are taken in the loop's grain only.
- `Clip: none` keeps the contour as given (warning if it leaves the grain): its cut crosses the GB.

Why clipping and skipping are needed: [Theory](theory.md#polycrystal).

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

- **Construction.** The first twin boundary (TB) is moved onto the nearest atomic plane; the atoms of
  the next m − 1 layers are moved in-plane onto the mirror lattice. No atom is added or removed.
  m is snapped to a multiple of 3 (minimum 3 planes) so that the second TB is also a coherent mirror
  ([Theory](theory.md#twin-lamellae)).
- **Extent.** A lamella spans its whole grain and ends at the GBs. In a single crystal it spans the
  periodic box, so the twin plane must be normal to a box axis (use `Orientation`, e.g. Z = [111]).
  Lamellae of one grain must be parallel and separated by at least one matrix layer.
- **Loops in lamellae.** A loop whose centre lies inside a lamella is read in the **twin frame**
  (`m`, `n`, `Side`, … are twin-crystal directions; R_twin is the parent turned by 180° about n) and
  snapped to the twin's planes. Every loop of a twinned grain or crystal is clipped at the TBs like at
  the GBs, and its matter slab only takes atoms of its own region.
- **BCC.** The pure-mirror {112} TB is unrelaxed; relax before MD.
- Output: an xyz column `twin` (0 = matrix, t = lamella t), `<name>_twins.txt` (plane, normal,
  TB point, m, thickness, atoms moved, R_twin), and a `region` column in `Poly_..._loops.txt`.

Examples: `Examples/fcc_twin_lamella.yaml` (single crystal, with loops in and around the lamellae),
`Examples/polycrystal_cu_twins.yaml`. For random lamellae in a polycrystal: `make_polycrystal_loops.py
--twins-per-grain` ([Input generators](generators.md)).

Checks with OVITO: [Validation](validation.md#generated-samples-ovito).

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

The grain assignment, interior crystallinity and orientations, GB distances, loops in grains and twins
were checked with OVITO (DXA, PTM, CNA); see [Validation](validation.md#generated-samples-ovito).
Building 2.2 M Cu atoms in 50 grains takes 0.2 s on 12 CPU threads.
