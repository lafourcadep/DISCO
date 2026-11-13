# Input generators

The scripts in `Examples/generators/` write DISCO YAML inputs for common configurations. Each
generated `<name>.yaml` comes with a **ground-truth sidecar** `<name>.json` that records what the sample
should contain. Compare it with the output of a dislocation analysis (DXA, DANCE, …) to validate that
analysis, or use it as the reference of a study.

| Script | Lattice | Writes |
|---|---|---|
| `make_single_loop.py` | BCC | one glide loop on any of the 48 slip systems: circle, square, rectangle, ⟨111⟩ parallelogram |
| `make_quadrupole.py` | BCC | four glide loops on two {110} planes: junction or annihilation; many in one sample (`--tiles`) |
| `make_prismatic_dipoles.py` | BCC | prismatic dipoles built from shear loops on a glide prism (rhombic or hexagonal), one or N |
| `make_prismatic_loops.py` | BCC | prismatic and mixed loops with atoms added or removed (interstitial, vacancy) |
| `make_fcc_loops.py` | FCC | perfect prismatic rhombi and Frank loops with atoms added or removed |
| `make_polycrystal_loops.py` | FCC, BCC | Voronoi polycrystal with a loop population and optional twin lamellae |
| `exec_input.py` | | shared helpers (no command line) |

## Common to all scripts

- Python 3 with numpy; run from any directory. `--help` lists every option.
- **Shared options**, with the same meaning everywhere: `--a0` (Å), `--format xyz|lammps|both|none`
  (the `Output_format` DISCO will use; default `lammps`), `-o/--output` (the sidecar goes next to it),
  and, for the random generators, `--seed` (default 0, so every run is reproducible).
- **Lengths** are in Å. `make_single_loop.py` and `make_quadrupole.py` also accept `--units a0`,
  which switches all their length options to units of a0.
- **Field method.** Every file uses `Method: line` and also writes a `Loopon_size`, so that switching
  to `Method: loopons` gives the grid the generator was designed for.
- **Mid-planes.** DISCO moves every loop centre onto a mid-plane between atomic planes
  ([Theory](theory.md#cut-placement)). The scripts reproduce this rule, record the snapped centre in the
  sidecar, and place their constructions (quadrupole faces, dipole faces) on mid-planes directly, so
  DISCO moves nothing.
- **Slip systems (BCC)** are read from `Parameters/BCC.txt` (`--bcc` to use another file): 48
  systems, 12 per ½⟨111⟩ Burgers vector, each listing 3 {110}, 3 {112} and 6 {123} planes. The
  DANCE analysis code uses the same file and numbering, so the sidecar's `expected_slip_system` can be
  compared directly with DANCE. FCC loops have no DANCE table: these fields are `null`.
- Defaults: BCC Ta (a0 = 3.3024 Å), FCC Cu (a0 = 3.615 Å).

## `make_single_loop.py`: one glide loop

One glide loop (**b** · **n** = 0) on a slip system of `BCC.txt` (`--system N`), on all systems of a
plane family (`--family 110|112|123|all`, one file each), or on any **m**, **n** (`--m … --n …`).

- `--shape circular` (radius `--radius`, default 8 a0), `square` (`--half`, 8 a0), `rectangular`
  (`--half-m` along **m**, screw sides, 16 a0; `--half-p` along **m** × **n**, edge sides, 8 a0;
  `--corner-radius`).
- `--shape parallelogram`: the faceted low-temperature BCC loop. `Dir_1` = **b** (screw sides,
  `--half-screw`, 16 a0) and `Dir_2` = the other ⟨111⟩ of the plane (mixed sides at 70.5° to **b**,
  `--half-mixed`, 8 a0). A {110} plane holds exactly two ⟨111⟩, so this shape exists on {110} only.
- Box: the smallest even cube leaving 12 a0 between the loop and the box faces (`--cells` to override).
  Default loopon step 0.4 a0. Files are named `circ_f112_sys04.yaml`, `rect_rounded_…`, `para_…`.

```bash
python3 make_single_loop.py --system 1                          # circle R = 8 a0, system 1
python3 make_single_loop.py --shape square --system 4           # {112}
python3 make_single_loop.py --shape rectangular --system 7 --corner-radius 13.2
python3 make_single_loop.py --shape parallelogram --family 110 --outdir db/   # the 12 {110} systems
python3 make_single_loop.py --family all --outdir db/           # all 48 systems
python3 make_single_loop.py --m 1 -1 1 --n 1 1 0 -o my_loop.yaml
```

Sidecar `expected`: one closed ½⟨111⟩ loop; for rectangles and parallelograms, the nominal screw, edge
or mixed line lengths (straight parts of the sides).

## `make_quadrupole.py`: junction and annihilation quadrupoles

Four glide loops on two {110} planes, whose corner reactions form a junction or annihilate; with
`--tiles`, many different quadrupoles in one sample. See [Quadrupoles](quadrupoles.md).

Sidecar `expected`: `closed_loops` is `null` (the loops meet on the corner lines and form a network);
`junctions` counts the corners expected to form an a⟨100⟩ junction; `corner_reactions` (or `tiles[]`)
gives, per corner, the two loops, the product Burgers vector, Frank's rule and the `BCC.txt` junction type.

## `make_prismatic_dipoles.py`: prismatic dipoles from shear loops

!!! warning "Not stable at finite temperature"
    These dipoles are built from shear loops only, so no atom is added or removed. The two loops of a
    dipole lie on one glide prism and attract each other: at finite temperature they glide together
    and annihilate, even without stress. For stable prismatic loops use `make_prismatic_loops.py`.

**Construction.** Shear loops with the same **b** sit on the faces of a glide prism around **b**:
two {110} planes (`--type rhombic`, 4 loops) or the three {110} planes containing **b** (`--type
hexagonal`, default, 6 loops). Every face normal points out of the tube, so the screw sides of
neighbouring faces coincide with opposite line senses and cancel. What remains are the edge sides at
the two tube ends: two prismatic loops of opposite sign, normal to **b**, separated by `--separation`
along **b**. The hexagonal ones have ⟨112⟩ sides, the glide prism of a BCC prismatic loop.

**Exact cancellation** needs the shared sides to coincide. The face planes are put on mid-planes
(this sets the side, `--side`), and the half-side and half-separation are written as
(k + ½) `Loopon_size`, with `Loopon_size` tuned close to `--loopon-size`. That is the exact contour
with `Method: line`, and the effective extent of the loopon grid with `Method: loopons`, so both methods
cancel the screw sides. The values used are printed.

**One dipole** (`--count 1`, default): centred at the origin, **b** = `--b` (default [1 −1 1]); for a
rhombic dipole, the planes `--n1 --n2` (default (110) and (011)). Box `--cells` (default 50).

**N dipoles** (`--count N`), drawn at random (`--seed`):

- **b**: one of the four ⟨111⟩ with a random sign (which end loop is positive). By default the four
  directions are dealt in shuffled rounds, so families stay balanced (12 dipoles = 3 per family);
  `--b-draw random` draws independently, `--b` fixes the direction, `--no-random-sign` the sign.
- Planes (rhombic): one of the three {110} pairs containing **b**. Separation: one value, or
  `--separation MIN MAX` for a uniform value per dipole.
- Position: a random BCC lattice site (cube corner or body centre), so every face stays on its
  mid-plane. A dipole is redrawn when its bounding sphere comes closer than `--gap` (10 Å) to another
  (periodic minimum image). With `--margin M`, every shear loop also stays M Å from the box faces, so no
  dipole crosses the periodic boundary (needed before opening a free surface, e.g. for indentation).
- Box: `--cells N` or `NX NY NZ` (unit cells), `--size LX LY LZ` (Å), or by default the cube in which
  the bounding spheres (with half the gap) fill `--fill` = 15 % of the volume.
- **One type and one side per sample**: the loopon step that cancels the screw sides is set by the side,
  and no step fits a rhombic and a hexagonal dipole at once (their half-sides are odd multiples of
  d/√3 and d/(2√3), d = {110} spacing). Only the separation varies between dipoles.

```bash
python3 make_prismatic_dipoles.py                                   # one hexagonal dipole
python3 make_prismatic_dipoles.py --type rhombic --side 40 --separation 120 --cells 60
python3 make_prismatic_dipoles.py --type rhombic --b 1 1 1 --n1 1 -1 0 --n2 0 1 -1
python3 make_prismatic_dipoles.py --count 10                        # 10 hexagonal dipoles, random b
python3 make_prismatic_dipoles.py --count 20 --type rhombic --side 40 --separation 60 120 --seed 3
# indentation slab: 12 dipoles (3 per <111> family), none crossing the box faces
python3 make_prismatic_dipoles.py --count 12 --side 40 --separation 50 --size 600 600 300 --margin 10 --gap 20
```

If the dipoles do not fit, the script stops and says how many it placed: enlarge the box, or lower
`--count`, `--fill`, `--gap` or `--margin`.

<div class="grid" markdown>
<figure markdown="span">
  ![Rhombic prismatic dipole](assets/prismatic.png){ height="220" }
  <figcaption>Rhombic prismatic dipole</figcaption>
</figure>
<figure markdown="span">
  ![Hexagonal prismatic dipole](assets/prismatic_hexagonal.png){ height="220" }
  <figcaption>Hexagonal prismatic dipole</figcaption>
</figure>
</div>

Sidecar `expected`: 2N closed ½⟨111⟩ loops, no screw part, and `per_dipole[]` (shear-loop numbers,
**b**, centre, planes, side, separation, centres of the two prismatic loops); `placement` records the
box, gap, margin, seed and **b** draw.

## `make_prismatic_loops.py`: BCC prismatic loops with matter

Prismatic loops that survive thermal MD: **b** is out of the loop plane (**b** · **n** ≠ 0) and DISCO
adds (interstitial) or removes (vacancy) the slab of atoms inside the loop
([Theory](theory.md#loops-with-b-n-0-adding-or-removing-matter)).

- `--habit 111`: pure-edge loop, **n** ∥ **b**, 3 {111} planes in the slab. `--habit 110`: mixed loop
  on a {110} plane at 35.3° from **b**, 1 plane.
- `--type interstitial | vacancy | mixed` (random per loop). `--b`: one ⟨111⟩ for all loops (default
  random per loop).
- `--shape circular | hexagonal`; `--radius` is the radius or the hexagon circumradius;
  `--corner-radius` rounds the hexagon. A regular hexagon's sides are parallel to its vertex
  directions. With `--habit 111` the first vertex is a ⟨112⟩ of the plane, so the **sides run along
  ⟨112⟩**, on the {110} glide prism of **b** (as in the dipoles and the polycrystal generator).
  `--vertex-dir` overrides it (it must lie in the plane of every loop, so give `--b` with it).
- With `--count N`, centres are random (`--seed`), bounding spheres at least `--gap` apart, and two loops
  of opposite type with the same **b** are never on overlapping glide prisms (they would glide together
  and annihilate). One loop without `--cells` is centred.
- `Loopon_size` defaults to 0.33 d (d: habit-plane spacing), the largest step that keeps the planes
  next to an interstitial cut smooth with `Method: loopons`.

```bash
python3 make_prismatic_loops.py                                     # 1 interstitial loop, R = 20 Å
python3 make_prismatic_loops.py --count 6 --type mixed --radius 15 --cells 50 --seed 2
python3 make_prismatic_loops.py --b 1 1 1 --type vacancy --habit 110 -o vac110.yaml
python3 make_prismatic_loops.py --shape hexagonal --radius 22 --corner-radius 2 -o hex.yaml
python3 make_prismatic_loops.py --shape hexagonal --b 1 1 1 --vertex-dir 1 -1 0 -o hex110.yaml  # <110> sides
```

Sidecar `expected`: one closed ½⟨111⟩ loop each, and `prismatic_loops[]` with type, normal and the
nominal atoms added (+) or removed (−), ρ A |**b** · n̂|.

## `make_fcc_loops.py`: FCC prismatic and Frank loops

FCC samples (`Lattice: FCC`, Cu by default) with loops made by adding or removing atoms: one loop, or
an irradiation-like population of interstitial and vacancy loops.

- `--kind prismatic` (default): perfect **b** = ½⟨110⟩, loop plane normal to **b** ({110}), 2 atomic
  planes. Exactly two {111} planes contain **b**, so the loop is a **rhombus** whose sides are its
  intersections with them: ⟨112⟩ directions at 70.5° / 109.5° (a DISCO `parallelogram`). Each side
  glides along **b** on its own {111} plane, so the loop moves rigidly along a parallelepiped glide
  prism: the FCC counterpart of the hexagonal ½⟨111⟩ prism in BCC. `--shape circular` gives a round loop.
- `--kind frank`: faulted **b** = a0/3⟨111⟩ on the {111} plane normal to **b**, sessile, one {111}
  layer added (extrinsic fault) or removed (intrinsic fault). `--shape hexagonal` (default) or
  `triangular`, with sides along ⟨110⟩ as in quenched and irradiated FCC metals, or `circular`.
- `--kind mixed`: each loop is Frank with probability `--frank-fraction` (default 0.5).
- `--type interstitial | vacancy | mixed`. `--size`: rhombus side, polygon circumradius or circle
  radius (Å). `--b`: one Burgers direction (a ⟨110⟩ for prismatic, a ⟨111⟩ for Frank loops).
- Placement as in `make_prismatic_loops.py` (perfect loops of opposite type and the same **b** never on
  overlapping glide prisms; Frank loops are sessile and only need `--gap`). Loops may cross the
  periodic boundary.
- `--element`, `--mass`, `--a0` for another metal (e.g. `--element Ni --a0 3.52`).

```bash
python3 make_fcc_loops.py                                          # 1 interstitial rhombus, side 25 Å
python3 make_fcc_loops.py --count 8 --type mixed --seed 2          # 8 prismatic rhombi
python3 make_fcc_loops.py --kind frank --shape triangular --type vacancy --size 20
python3 make_fcc_loops.py --kind mixed --count 12 --type mixed --frank-fraction 0.5 --cells 60
python3 make_fcc_loops.py --kind prismatic --b 1 -1 0 --size 30 --element Ni --a0 3.52
```

Sidecar `expected.loops[]`, per loop: kind, type, Burgers vector (`1/2[1 1 0]`, `1/3[1 1 1]`), number of
atomic planes, nominal atoms added (ρ A |**b** · n̂|), and the glide planes (prismatic) or the stacking
fault (Frank).

## `make_polycrystal_loops.py`: dislocation-populated polycrystal

Draws the grains (seeds uniform in the box, uniformly random orientations), writes them as an explicit
`Grain_list` (Center + Bunge Euler), then adds loops until the nominal line length reaches `--density`
(m⁻²) or `--loops`. See [Polycrystals](polycrystal.md) and
[Theory](theory.md#polycrystal) for clipping and skipping.

- Box `--box L` or `LX LY LZ` (Å), `--grains`, `--lattice FCC|BCC` (Cu or Ta defaults).
- Each loop: centre uniform in the box, grain = nearest seed (grains get loops in proportion to their
  volume), bounding spheres at least `--gap` apart (periodic).
- `--kinds name:weight,...`:

    | Lattice | Kind | Loop |
    |---|---|---|
    | FCC | `glide` | ½⟨110⟩{111}, circular |
    | FCC | `shockley` | ⅙⟨112⟩{111}, circular, stable sense only |
    | FCC | `frank` | ⅓⟨111⟩, hexagon with ⟨110⟩ sides |
    | FCC | `prismatic` | ½⟨110⟩ rhombus on its {111} glide prism |
    | BCC | `glide` | ½⟨111⟩{110}, circular |
    | BCC | `prismatic` | ½⟨111⟩ hexagon with ⟨112⟩ sides on its {110} glide prism |

    Default FCC `glide:4,frank:3,prismatic:2,shockley:1`, BCC `glide:3,prismatic:2`.
- `--size` (median radius or reach, Å) and `--size-sigma` (log-normal log-sd); `--type` for the loops with
  **b** · **n** ≠ 0.
- **Twins**: `--twins-per-grain λ` gives each grain a Poisson(λ) number of parallel lamellae (random
  {111} FCC / {112} BCC variant), thickness log-normal (`--twin-thickness`, `--twin-sigma`), separated
  by at least 3 planes plus up to `--twin-gap` Å; a stack stops where it would leave its grain. Loops
  whose centre falls in a lamella are read by DISCO in the twin frame. `--loops 0` writes twins only.

```bash
python3 make_polycrystal_loops.py                                  # Cu, 8 grains, 1e16 m^-2
python3 make_polycrystal_loops.py --grains 12 --box 200 --density 5e16 --seed 5
python3 make_polycrystal_loops.py --kinds glide:1,frank:1 --size 25 --type interstitial
python3 make_polycrystal_loops.py --lattice BCC --kinds glide:2,prismatic:1
python3 make_polycrystal_loops.py --twins-per-grain 1.5 --twin-thickness 15 --density 2e16
python3 make_polycrystal_loops.py --twins-per-grain 2 --loops 0    # twins only
```

The sidecar lists the grains (centre, Euler, R), the twins (grain, plane, point, R_twin) and the loops
(grain, kind, type, crystal and box Burgers vectors, approximate region, `may_clip` when the bounding
sphere leaves the grain cell). Its `expected` block counts the loops *before* clipping; DISCO's
`Poly_..._loops.txt` gives the final status of each loop and the density after clipping.

## Sidecar (`<name>.json`), schema 1

| Key | Content |
|---|---|
| `schema` | sidecar format version (1) |
| `generator`, `command` | the script and the exact command line that wrote the file |
| `general` | the YAML `general:` block (polycrystal: lattice, element, a0, box) |
| `loops[]` | one record per input loop (below) |
| `expected` | what the sample should contain (below) |
| script-specific | `quadrupole`, `tiling`, `placement`, `grains`, `twins`, `note` |

**Per-loop record** (all scripts but the polycrystal one): shape and sizes; `center_requested` and
`center_snapped` (after the mid-plane rule); `m`, `n`, `burgers_A`, `plane_family`;
`expected_burgers_type` and `expected_slip_system` (indices in `BCC.txt`, `null` for FCC); `loopons`,
the number of loopons DISCO would generate with `Method: loopons`.

**`expected`**, keys present in every sidecar:

| Key | Meaning |
|---|---|
| `closed_loops` | number of closed loops a dislocation analysis should find; `null` when the loops form a network (quadrupoles) |
| `burgers_families` | Burgers vector families present, e.g. `["1/2<111>"]`, `["1/2<110>", "1/3<111>"]` |
| `junctions` | number of junctions expected (quadrupoles: Frank-favourable a⟨100⟩ corners) |

plus, per script: nominal line lengths (single loop), `corner_reactions` / `tiles` (quadrupoles),
`per_dipole` (dipoles), `prismatic_loops` (BCC prismatic loops), `loops` (FCC loops),
`nominal_density_m2` (polycrystal).

Coordinates follow the input file: origin at the box centre ([Conventions](conventions.md#origins)).
