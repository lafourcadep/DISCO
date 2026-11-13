# Conventions

This page fixes the notation and the sign conventions used everywhere else: in the input file, in
the code, in the generators and in their sidecars. The other pages link here instead of restating them.

## Notation

| Symbol | Meaning |
|---|---|
| a0 | cubic lattice parameter (Å) |
| **m** | Burgers vector in Miller-type components, as written in the input (`m:`): integers, or exact fractions p/q for FCC partials |
| **b** = (a0/2) **m** | Burgers vector (Å); \|b\| = a0 \|m\| / 2 |
| m̂ = **b** / \|b\| | unit Burgers direction (used in the formulas of [Theory](theory.md)) |
| **n**, n̂ | loop-plane normal: integer Miller indices (`n:`), and the unit vector along it |
| d | spacing of the atomic planes of the family **n** ([Theory](theory.md#cut-placement)) |
| ν | Poisson's ratio, fixed at 0.35 (`KROUPA_NU` in `Code/kroupa.hpp`) |
| ρ | atomic density: 2/a0³ (BCC), 4/a0³ (FCC) |
| A | area of the loop contour |
| R | rotation of a crystal (grain, oriented box, twin lamella): **x**_box = R **x**_crystal |

Families are written in the usual way: ½⟨111⟩ means **b** = (a0/2)⟨111⟩ (**m** = ⟨111⟩),
⅓⟨111⟩ means **m** = ⅔⟨111⟩, ⅙⟨112⟩ means **m** = ⅓⟨112⟩.

## Loop orientation and displacement jump

- The loop contour C runs **counter-clockwise about n** (right-handed). This line sense is the
  one of the dislocation.
- The cut S is the flat surface bounded by C. Across it, the displacement jumps by
  **u⁺ − u⁻ = −b**, where + is the side n̂ points to.
- Reversing **n** (or **m**) reverses the sign of the loop.

## Loop character

| Condition | Loop | `Matter: auto` |
|---|---|---|
| **m** · **n** = 0 | glide (shear) loop: **b** in the loop plane | no atom added or removed |
| **m** · **n** > 0 | vacancy loop (prismatic if **m** ∥ **n**, mixed otherwise) | slab of thickness \|b · n̂\| removed |
| **m** · **n** < 0 | interstitial loop | slab of thickness \|b · n̂\| duplicated |

For FCC faulted loops:

- **Frank** (**m** = ±⅔⟨111⟩ ∥ **n**): vacancy (m · n > 0) leaves an intrinsic fault, interstitial
  (m · n < 0) an extrinsic fault.
- **Shockley** (**m** = ⅓⟨112⟩ in the {111} plane): only one sense of **b** makes a stable
  intrinsic fault. On each {111} plane three of the six ⅓⟨112⟩ are stable; for
  **n** = [111], **m** = [−⅓, −⅓, ⅔] is stable and [⅓, ⅓, −⅔] is not. DISCO warns for the wrong sense.

## Frames

- **Crystal frame**: `m`, `n`, `Side`, `Dir_1`, `Dir_2`, `Vertex_dir` and the `Vertices` offsets
  are crystal directions: of the cubic crystal, of the oriented crystal (`Orientation`), of the
  loop's grain (`Grain: k`), or of the twin lamella that contains the loop centre.
- **Box (lab) frame**: `Center`, the box, the twin-plane `Center`, and all output positions.
- **Rotation R**: its rows are the crystal directions of the box axes x, y, z, so
  **x**_box = R **x**_crystal. `Orientation: {X, Y, Z}` gives these rows directly
  (normalised); `Euler: [φ1, Φ, φ2]` gives them through the Bunge convention, g = Rᵀ.
- **Twin frame**: R_twin = R · (2 n̂n̂ᵀ − I), the parent turned by 180° about the twin-plane
  normal (Σ3: equivalent to 60° about ⟨111⟩).

## Origins

- **Input** (loop centres, twin points, grain seeds) and the generators' sidecars: origin at the
  **box centre**.
- **Output** files (`.xyz`, `.lmp`): origin at a **box corner**; every position is shifted by
  +L/2 (half of each box length). Add this shift to compare a sidecar position with a DXA result.

## Units

Every length is in Å: a0, `Center`, sizes, `Loopon_size`, `Segment_tolerance`, displacements.
**n** is dimensionless (normalised internally). Masses are in Dalton. Angles are in degrees.
