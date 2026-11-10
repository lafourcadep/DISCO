# Input file

## Running

```bash
./exec_cuda <input.yaml>          # GPU
./exec_cuda --cpu <input.yaml>    # CPU backend of the same binary
OMP_NUM_THREADS=16 ./exec_cpu <input.yaml>
```

`exec_cuda` falls back to the CPU backend, with a warning, when it finds no usable CUDA device.
The first lines of the log name the backend (`--- GPU: <name> ---` or
`--- CPU: <N> OpenMP threads ---`).

At the end the program prints a timing summary. It gives the wall-clock time and the share of each phase
(input, lattice, loopons, transfers, kernels, final positions, output), then the kernel time
per loop and its throughput in atom–loopon pairs per second. The per-loop lines appear only when
there are 20 loops or fewer.

## Format

The input is YAML: a `general:` block and a `dislocations:` list with one entry per loop.
Key order is free and `#` comments are allowed. **An unknown or missing key stops the program
with a message**, so a typo cannot silently fall back to a default.

!!! note "Units"
    Every length is in Ångström.

```yaml
general:
  a0: <lattice_parameter>                # Å
  Nx: <cells_x>                          # box = Nx × Ny × Nz cubic cells
  Ny: <cells_y>
  Nz: <cells_z>
  Loopon_size: <step>                    # Å, loopon grid step, the same for every loop
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
    m: [mx, my, mz]          # Burgers vector direction (integer Miller indices)
    n: [nx, ny, nz]          # loop-plane normal (integer Miller indices)
    Matter: auto             # optional, loops with b·n ≠ 0 only: auto (default) adds or removes
                             #   atoms, none applies the displacement field only
  - ...
```

The Burgers vector is **b** = (a0/2) **m**, so |b| = a0 |m| / 2. The log prints the size of each loop
in Å and in units of |b|.

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

When m is not in the loop plane, the loop is prismatic (m ∥ n) or mixed. Opening or closing the cut
by |b · n̂| needs matter. With `Matter: auto`, the code adds or removes the atoms of a slab of
|b · n̂| = k d inside the contour. k is an integer because b is a lattice vector.

| sign of m · n | loop | atoms |
|---|---|---|
| m · n > 0 | vacancy | slab atoms removed before the displacement sum |
| m · n < 0 | interstitial | slab atoms duplicated; each copy is displaced by u(x) + sign(s) b and lands on the other lip of the cut |

This convention follows the Kroupa field: across the cut, u⁺ − u⁻ = −b.
The log prints the number of atoms added or removed and the expected value ρ A |b · n̂|.
The box is not deformed. The volume change ±A |b · n̂| is printed, so relax the box with NPT or
`fix box/relax`.

- The slab is bounded by mid-planes on either side of the cut. It holds 3 {111} planes for
  ½⟨111⟩ on (111), 1 {110} plane for ½⟨111⟩ on {110}, and 2 {100} planes for a⟨100⟩ on (100).
- For an interstitial loop, the atoms kept next to the cut are only d/2 from the loopon plane.
  Use `Loopon_size` ≤ 0.35 d, for example 0.33 Å for {111}. The code warns otherwise.
- For m ∥ n, the in-plane grid axes are n × e_k (e_k is the Cartesian axis with the smallest |n_k|)
  and that vector × n. `Half_m` of a rectangle is measured along the first one, and the first vertex
  of a regular polygon without `Vertex_dir` lies along it.
- Any shape works, since the slab is cut by the loop contour. Sharp polygon corners need no rounding:
  for interstitial and vacancy hexagons on (111) (R = 22 Å, `Loopon_size` 0.33 Å), adaptive CNA finds
  no extra non-BCC atoms at the corners, and DXA finds one closed ½⟨111⟩ loop, with or without
  `Corner_radius`.
- `Examples/generators/make_prismatic_loops.py` writes such inputs.

A loop made this way is stable at temperature without stress. It can glide along b but cannot
vanish without climb. In contrast, the dipoles built from shear loops (`make_prismatic_dipole.py`)
contain no extra matter. Their two loops share one glide prism and annihilate by glide.

## Loop centre and cut plane

The loop surface must lie halfway between two dense atomic planes of its (hkl) = n family.
The code moves the requested centre along n onto the nearest such mid-plane, without changing the
in-plane position. It prints the shift and the plane spacing d:

- d = a0/|n| if h + k + l is even;
- d = a0/(2|n|) if h + k + l is odd.

This gives d = 2.335 Å for {110}, 1.348 Å for {112} and 0.883 Å for {123} (a0 = 3.3024 Å).
A centre exactly on an atomic plane goes to the mid-plane on the +n side.

## Example: one circular and one rounded-rectangle loop

```yaml
general:
  a0: 3.3024
  Nx: 60
  Ny: 60
  Nz: 60
  Loopon_size: 1.32096
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

More templates are in `Examples/` (`cubic_base.yaml`, `orthorhombic_base.yaml`,
`quad_junc.yaml`, `quad_anni.yaml`, `quadrupole_circle.yaml`, `rectangular.yaml`).
The [generators](generators.md) write such files for you.

## Units

| Quantity | Unit |
|---|---|
| Lattice parameter `a0` | Å |
| Atom and loopon positions, all input lengths | Å |
| Burgers vector b | Å |
| Loop-plane normal n | dimensionless (normalised internally) |
| Loopon area dS | Å² |
