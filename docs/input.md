# Input file

## Running

```bash
./exec_cuda <input.yaml>          # GPU
./exec_cuda --cpu <input.yaml>    # CPU backend of the same binary
OMP_NUM_THREADS=16 ./exec_cpu <input.yaml>
./exec_mvcc --fp32 <input.yaml>   # Apple GPU (see CPU or GPU)
```

`exec_cuda` falls back to the CPU backend, with a warning, when it finds no usable CUDA device.
The first lines of the log name the backend (`--- GPU: <name> ---` or
`--- CPU: <N> OpenMP threads ---`).

At the end the program prints a timing summary. It gives the wall-clock time and the share of each phase
(input, lattice, loopons, transfers, kernels, final positions, output), then the kernel time
per loop and its throughput in atom–segment (`Method: line`) or atom–loopon (`Method: loopons`)
pairs per second. The per-loop lines appear only when
there are 20 loops or fewer.

## Format

The input is YAML: a `general:` block and a `dislocations:` list with one entry per loop.
Key order is free and `#` comments are allowed. **An unknown or missing key stops the program
with a message**, so a typo cannot silently fall back to a default.

!!! note "Units"
    Every length is in Ångström.

```yaml
general:
  Lattice: <BCC|FCC>                     # optional, default BCC
  Element: <symbol>                      # optional, default Ta: species name in the output files
  Mass: <mass>                           # optional (Dalton), LAMMPS Masses; default: standard mass
                                         #   of Element (common BCC / FCC metals), required otherwise
                                         #   A warning is printed if Element is not of this Lattice
                                         #   at room temperature (e.g. Ni with BCC), or if Lattice is
                                         #   FCC and Element is left at its default (Ta)
  a0: <lattice_parameter>                # Å
  Orientation: {X: [..], Y: [..], Z: [..]}  # optional: crystal directions of the box axes
                                         #   (see "Oriented crystal"); default [100] [010] [001]
  Nx: <cells_x>                          # box = Nx × Ny × Nz cubic cells, or with Orientation
                                         #   Nx × Ny × Nz lattice periods along X, Y, Z
  Ny: <cells_y>
  Nz: <cells_z>
  Method: <line|loopons>                 # optional, default line (see "Field method")
  Segment_tolerance: <eps>               # line only, optional, default 0.01 Å
  Loopon_size: <step>                    # Å, loopon grid step, the same for every loop;
                                         #   required with loopons, optional (unused) with line
  Output_format: <xyz|lammps|both|none>  # optional, default xyz; none = no file (timing runs)

dislocations:
  - Shape: <circular|rectangular|parallelogram|polygon|hexagonal>
    Radius: <value>          # circular only (Å)
    Half_m: <value>          # rectangular only (Å): half-length along m (screw sides)
    Half_p: <value>          # rectangular only (Å): half-length along m × n (edge sides)
    Corner_radius: <value>   # rectangular (required), polygon / hexagonal (optional, default 0), Å;
                             #   0 = sharp corners
    Side: [sx, sy, sz]       # rectangular only, optional: direction of the Half_m axis
                             #   (in the plane, Side·n = 0); default m. Half_p is along Side × n
    Half_1: <value>          # parallelogram only (Å): half-length along Dir_1
    Half_2: <value>          # parallelogram only (Å): half-length along Dir_2
    Dir_1: [d1x, d1y, d1z]   # parallelogram only: first side direction, in the plane
    Dir_2: [d2x, d2y, d2z]   # parallelogram only: second side direction, in the plane, not ∥ Dir_1
    Sides: <N>               # polygon only: regular polygon with N >= 3 sides (hexagonal: 6)
    Circumradius: <value>    # polygon / hexagonal (Å): centre-to-vertex distance
    Vertex_dir: [vx, vy, vz] # polygon / hexagonal, optional: centre -> first vertex, in the plane;
                             #   default: the in-plane grid axis (see below)
    Vertices: [[x, y, z], ...]  # polygon only, instead of Sides / Circumradius / Vertex_dir:
                             #   convex contour, offsets from Center in Å, in the plane, in order
    Center: [cx, cy, cz]     # Å, origin at the box centre (moved onto a mid-plane, see below)
    m: [mx, my, mz]          # Burgers vector b = (a0/2) m: integers, or exact fractions p/q
                             #   for FCC partials (Frank [2/3, 2/3, 2/3], Shockley [1/3, 1/3, -2/3])
    n: [nx, ny, nz]          # loop-plane normal (integer Miller indices)
    Matter: auto             # optional, loops with b·n ≠ 0 only: auto (default) adds or removes
                             #   atoms, none applies the displacement field only
  - ...
```

The Burgers vector is **b** = (a0/2) **m**, so |b| = a0 |m| / 2. The log prints the size of each loop
in Å and in units of |b|. Signs, frames and origins are defined in [Conventions](conventions.md).

## Lattice

| `Lattice` | Atoms per cubic cell | Perfect Burgers vectors (b = a0/2 m) | Typical glide plane |
|---|---|---|---|
| `BCC` (default) | 2 | m = ⟨111⟩ (½⟨111⟩), m = ⟨200⟩ (⟨100⟩) | {110}, {112}, {123} |
| `FCC` | 4 | m = ⟨110⟩ (½⟨110⟩), m = ⟨200⟩ (⟨100⟩) | {111} |

The program warns when b = (a0/2) m is not a translation of the lattice (for example m = ⟨111⟩ in
FCC): the cut surface of such a loop leaves a fault and the atoms added or removed by a prismatic
loop would not restore the perfect crystal. The two faulted FCC loops below are recognised and get
an information line instead. Perfect FCC loops dissociate into Shockley partials when the sample is
relaxed (energy minimization or MD), as expected.

### Faulted loops (FCC)

Write the fractional m as exact fractions (`2/3`, `-1/3`, …); a decimal such as `0.6667` stops the
program.

| Loop | m | n | Inside the contour |
|---|---|---|---|
| Frank, vacancy | ⅔⟨111⟩, m·n > 0 | ∥ m ({111}) | one {111} layer removed: intrinsic stacking fault |
| Frank, interstitial | ⅔⟨111⟩, m·n < 0 | ∥ m ({111}) | one {111} layer inserted: extrinsic stacking fault |
| Shockley | ⅓⟨112⟩, m·n = 0 | {111} containing m | glide loop: intrinsic stacking fault |

Inserted Frank layers are shifted onto the third stacking position, and only one sense of a Shockley
**b** gives a stable fault: for n = [1 1 1], m = [-1/3, -1/3, 2/3] is stable and m = [1/3, 1/3, -2/3]
is not (DISCO warns; reverse m or n). See [Theory](theory.md#faulted-loops-fcc).

Example: `Examples/fcc_faulted_loops.yaml` (Cu: interstitial and vacancy Frank loops, one Shockley
loop). DXA (FCC) finds ⅓⟨111⟩ Frank and ⅙⟨112⟩ Shockley loops, with HCP (fault) layers inside.

Example (Cu, two glide loops on {111}):

```yaml
general:
  Lattice: FCC
  Element: Cu                # Mass 63.546 from the built-in table
  a0: 3.615
  Nx: 40
  Ny: 40
  Nz: 40
dislocations:
  - Shape: rectangular
    Half_m: 30.0
    Half_p: 22.0
    Corner_radius: 0.4
    Center: [-25, -25, 0]
    m: [1, -1, 0]
    n: [1, 1, 1]
  - Shape: hexagonal
    Circumradius: 20.0
    Corner_radius: 2.0
    Center: [25, 25, 10]
    m: [1, 0, 1]
    n: [1, 1, -1]
```

`Examples/generators/make_fcc_loops.py` writes FCC inputs (prismatic and Frank loops); see
[Input generators](generators.md).

### Polycrystal

A top-level `polycrystal:` block (box, number of grains or explicit grains) replaces `Nx`, `Ny`, `Nz`
and `Orientation`; see [Polycrystals](polycrystal.md).

### Twins

A top-level `twins:` list adds coherent twin lamellae ({111} FCC, {112} BCC) to an oriented single
crystal (twin plane normal to a box axis) or to polycrystal grains; loops inside a lamella use the
twin frame. See [Polycrystals: Twin lamellae](polycrystal.md#twin-lamellae).

### Oriented crystal

`Orientation` gives the crystal directions of the box axes, for example the usual FCC dislocation
setup:

```yaml
general:
  Lattice: FCC
  Element: Cu
  a0: 3.615
  Orientation: {X: [1, -1, 0], Y: [1, 1, -2], Z: [1, 1, 1]}   # (111) planes horizontal
  Nx: 32            # periods along X: 32 × a0/√2
  Ny: 18            # periods along Y: 18 × a0 √6 / 2
  Nz: 13            # periods along Z: 13 × a0 √3
```

- X, Y, Z are integer directions, orthogonal and right-handed (Z along X × Y); otherwise the program
  stops.
- `Nx`, `Ny`, `Nz` count **lattice periods** along X, Y, Z, the shortest lattice translation along each
  direction ([Theory](theory.md#oriented-crystal)), so the box is exactly periodic. The log prints
  each period.
- **Frames.** Crystallographic vectors stay in the crystal frame: `m`, `n`, `Side`, `Dir_1`, `Dir_2`,
  `Vertex_dir` and the `Vertices` offsets. `Center` is in the lab frame (box axes, origin at the box
  centre). The default in-plane axis of a shape (no `Side` / `Vertex_dir`) is chosen in the crystal
  frame, so a shape keeps its crystallographic orientation in any box.
- Everything else (mid-plane snap, matter slab, Frank and Shockley rules, field) works as in the
  cube-aligned box. Without `Orientation` the output is unchanged.

Example: `Examples/fcc_oriented_loops.yaml`. Checks: [Validation](validation.md#generated-samples-ovito).

## Field method

The displacement field of each loop is computed in one of two ways
(see [Theory](theory.md#line-integral-method-line)).

| `Method` | Field | Resolution parameter |
|---|---|---|
| `line` (default) | Closed-form line integral: exact isotropic field of the contour, made of straight segments | `Segment_tolerance` ε: circles and rounded corners are cut into chords of sagitta ≤ ε |
| `loopons` | Kroupa sum over a grid of infinitesimal loops tiling the cut | `Loopon_size`: grid step |

- Straight sides (rectangles, parallelograms, polygons) are exact with `line`; circles and rounded
  corners are cut into chords of sagitta ≤ ε (a circle's chords keep its area; a rounded corner loses
  at most ~ε × arc length).
- `loopons` approximates the contour by the grid; the error falls with `Loopon_size` and is largest at
  the cores. `line` is exact and, having far fewer segments than a loop has loopons, faster
  ([Theory](theory.md#loopon-sum-method-loopons), [CPU or GPU](backends.md#speed)).
- `line` became the default after validation. Inputs written before it, with only `Loopon_size`,
  now run with `line`; add `Method: loopons` to reproduce older samples exactly.

## Loop shapes

**Circle.** A disc of radius `Radius` in the plane normal to n.

**Rectangle.** By default the sides run along m (screw character, half-length `Half_m`) and along
m × n (edge character, half-length `Half_p`). `Corner_radius` rounds the corners.
`Side` turns the rectangle in its plane without changing b. `Half_m` is then measured along
`Side` and `Half_p` along `Side × n`, and the sides have mixed character. The quadrupole generator
uses this to put the sides of neighbouring loops exactly on the line where their planes meet.

**Parallelogram.** Sides along `Dir_1` and `Dir_2`, which need not be orthogonal, so the
characters of the two side pairs are set independently. The loopon grid is oblique, with step
`Loopon_size` along each direction, and each loopon carries the cell area
`Loopon_size² · sin(Dir_1, Dir_2)`. The tiles follow both sides exactly. The effective half-lengths
are (k + ½) `Loopon_size` with k = floor(Half / `Loopon_size`).

A typical use is the faceted shape of a low-temperature BCC loop on {110}: `Dir_1` = b (screw sides),
`Dir_2` = the other ⟨111⟩ of the plane (mixed sides at 70.5° to b). For example, for b = ½[111] on (1−10):
`Dir_1: [1, 1, 1]`, `Dir_2: [1, 1, -1]`
(`Examples/generators/make_single_loop.py --shape parallelogram` writes these).

<figure markdown="span">
  ![Parallelogram loop](assets/parallelogram.png){ width="300" }
  <figcaption>Faceted ⟨111⟩ loop: screw sides along b, mixed sides at 70.5°.</figcaption>
</figure>

**Polygon / hexagon.** `Shape: hexagonal` is a regular hexagon. `Shape: polygon` is a regular
polygon with `Sides` sides, or any strictly convex polygon given by `Vertices`. A regular polygon is
set by `Circumradius` (centre to vertex) and `Vertex_dir`, the in-plane direction of its first vertex.
Without `Vertex_dir`, the first vertex lies along the default grid axis: m projected on the plane or,
for m ∥ n, n × e_k (e_k is the Cartesian axis with the smallest |n_k|). For n = ⟨111⟩ this axis is a
⟨110⟩ direction, so the sides run along ⟨112⟩. Give `Vertex_dir` along a ⟨112⟩ of the plane to get
sides along ⟨110⟩.
`Corner_radius` rounds every corner with an arc tangent to both sides. It is clamped to the largest
value that keeps every side, which for a regular polygon is the apothem; at that value the polygon
becomes its inscribed circle. The loopon grid is square, aligned with `Vertex_dir` (or the default
axis), and nodes on the contour are included. The log prints the exact contour area next to the
loopon area. `Vertices` must be strictly convex (no three aligned) and listed in order around the
contour, in either direction.

A typical use is the ½⟨111⟩{111} prismatic loop of BCC metals (interstitial or vacancy), for example:

```yaml
  - Shape: hexagonal
    Circumradius: 22.0
    Vertex_dir: [1, -1, 0]   # sides along <112>
    Corner_radius: 2.0       # optional
    Center: [0, 0, 0]
    m: [-1, -1, -1]          # m . n < 0: interstitial
    n: [1, 1, 1]
```

## Prismatic and mixed loops (b · n ≠ 0)

When m is not in the loop plane, the loop is prismatic (m ∥ n) or mixed. With `Matter: auto` (default)
DISCO adds or removes the atoms of the slab of thickness |b · n̂| inside the contour, so the loop
survives thermal MD ([Theory](theory.md#loops-with-b-n-0-adding-or-removing-matter)).
`Matter: none` applies the displacement field only.

| sign of m · n | loop | atoms |
|---|---|---|
| m · n > 0 | vacancy | slab atoms removed before the displacement sum |
| m · n < 0 | interstitial | slab atoms duplicated; each copy is displaced by u(x) + sign(s) b and lands on the other lip of the cut |

- The log prints the number of atoms added or removed and the expected value ρ A |b · n̂|, and the
  volume change ±A |b · n̂|. The box is not deformed: relax it with NPT or `fix box/relax`.
- With `Method: loopons`, use `Loopon_size` ≤ 0.35 d (e.g. 0.33 Å for BCC {111}); the code warns
  otherwise. `Method: line` has no such limit.
- For m ∥ n, the in-plane grid axes are n × e_k (e_k is the Cartesian axis with the smallest |n_k|)
  and that vector × n. `Half_m` of a rectangle is measured along the first one, and the first vertex
  of a regular polygon without `Vertex_dir` lies along it.
- Any shape works, and sharp polygon corners need no rounding.
- `Examples/generators/make_prismatic_loops.py` (BCC) and `make_fcc_loops.py` (FCC) write such inputs.

## Loop centre and cut plane

The loop surface must lie halfway between two dense atomic planes of its (hkl) = n family. DISCO
moves the requested centre along n onto the nearest such mid-plane, without changing its in-plane
position, and prints the shift and the plane spacing d (rule and values in
[Theory](theory.md#cut-placement)). A centre exactly on an atomic plane goes to the mid-plane on the
+n side.

## Example: one circular and one rounded-rectangle loop

```yaml
general:
  a0: 3.3024
  Nx: 60
  Ny: 60
  Nz: 60
  Method: line
  Output_format: lammps

dislocations:
  - Shape: circular
    Radius: 26.4192
    Center: [-39.6288, -39.6288, 0]
    m: [1, -1, 1]
    n: [1, 1, 0]
  - Shape: rectangular
    Half_m: 34.3195
    Half_p: 19.8144
    Corner_radius: 9.9072
    Center: [39.6288, 39.6288, 0]
    m: [-1, 1, 1]
    n: [1, 0, 1]
```

More templates are in `Examples/`: BCC (`cubic_base.yaml`, `orthorhombic_base.yaml`,
`quad_junc.yaml`, `quad_anni.yaml`, `quadrupole_circle.yaml`, `rectangular.yaml`), FCC
(`fcc_loops.yaml`, `fcc_faulted_loops.yaml`, `fcc_oriented_loops.yaml`, `fcc_twin_lamella.yaml`) and
polycrystals (`polycrystal_cu.yaml`, `polycrystal_cu_loops.yaml`, `polycrystal_cu_twins.yaml`).
The [generators](generators.md) write such files for you.

## Units

Every length is in Å, including `Segment_tolerance` and `Loopon_size`; n is dimensionless and
normalised internally ([Conventions](conventions.md#units)).
