# Polycrystals

A `polycrystal:` block replaces `Nx`, `Ny`, `Nz` (and `Orientation`) and builds a periodic
Voronoi polycrystal: each grain has a seed and an orientation, and each atom belongs to the grain
of its nearest seed (periodic minimum image).

!!! note "Status"
    Grains and grain boundaries (GBs) are built. Dislocation loops inside grains are the next step:
    for now a `polycrystal` input has no `dislocations` list.

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

## Output

- `Poly_<G>g_<Lx>x<Ly>x<Lz>.xyz` / `.lmp` (box lengths rounded to Å). The extended XYZ file has an
  extra per-atom column `grain` (1-based), which OVITO reads as a particle property. The LAMMPS file
  keeps `atom_style atomic` (one atom type); use the xyz file, or the table below, for grain IDs.
- `Poly_..._grains.txt`: one line per grain: id, seed (box-centred, Å), Bunge Euler angles (degrees),
  the 9 entries of R (rows: crystal directions of the lab x, y, z axes; lab = R · crystal) and the
  number of atoms.
- The log prints the same table.

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
