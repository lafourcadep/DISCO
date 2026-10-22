#!/usr/bin/env python3
"""
replicate_quadrupole.py
=======================
Scale, replicate, and optionally randomise the quadrupole of dislocation loops
defined in quadrupole_circle.yaml, producing a ready-to-run input file for the
Kroupa CUDA code (exec_cuda).

Base quadrupole (quadrupole_circle.yaml)
---------------------------------------
  Box:    40 × 40 × 40 unit cells (= 40 × 40 × 40 a0 in internal coordinates)
  Radius: 8 a0
  4 loops at centres (±6, ±6, 0) a0 — COM at origin.

  Loop 1  m = [ 1,-1,-1]  n = [ 1, 1, 0]  centre (-6,-6, 0)
  Loop 2  m = [-1,-1, 1]  n = [-1, 1, 0]  centre (-6,+6, 0)
  Loop 3  m = [-1,-1,-1]  n = [-1, 1, 0]  centre (+6,-6, 0)
  Loop 4  m = [ 1,-1, 1]  n = [ 1, 1, 0]  centre (+6,+6, 0)

What this script does
---------------------
  1. Scale  — change the loop radius; the centre offsets from the COM scale
               proportionally so the quadrupole geometry is preserved.
  2. Replicate — tile the quadrupole Nrx × Nry × Nrz times along x, y, z.
                 The simulation box grows accordingly.
  3. Noise  — add an independent random displacement to the COM of each tile,
               simulating positional disorder of the quadrupoles.
               Each set of 4 loops within a tile is shifted rigidly.

Usage
-----
  python replicate_quadrupole.py [options]

  All lengths on the command line are in Angstrom (a0 = 3.3024 Å by default).

  # Keep original radius, 2×2×1 tiling, no noise
  python replicate_quadrupole.py --nrx 2 --nry 2 -o quad_2x2x1.yaml

  # Scale radius to 12 a0 (39.6288 Å), 3×3×3 tiling
  python replicate_quadrupole.py --radius 39.6288 --nrx 3 --nry 3 --nrz 3 -o quad_r12_3x3x3.yaml

  # With Gaussian noise (std = 2 a0 = 6.6048 Å) on the COM of each tile, in xy only
  python replicate_quadrupole.py --nrx 4 --nry 4 --noise-amp 6.6048 --noise-axes xy --seed 42

  # Full example
  python replicate_quadrupole.py --radius 33.024 --nrx 2 --nry 2 --nrz 1 \\
      --noise-amp 4.9536 --noise-type gaussian --noise-axes xyz --seed 0 \\
      --format both -o quad_custom.yaml
"""

import argparse
import math
import random
import sys

# =============================================================================
# Base quadrupole template  (from quadrupole_circle.yaml)
# All positions are in units of a0.  Origin = box centre.
# =============================================================================

BASE_A0           = 3.3024   # Å
BASE_CELL_SIZE_X  = 50       # unit cells per axis  (= 40 a0)
BASE_CELL_SIZE_Y  = 40       # unit cells per axis  (= 40 a0)
BASE_CELL_SIZE_Z  = 32       # unit cells per axis  (= 40 a0)
BASE_RADIUS       = 8.0      # a0
BASE_LOOPON_SIZE  = 0.4      # a0 (same grid as the base file)

# Centre offsets relative to the tile COM (which sits at the origin for the
# single-tile case).  Scale factor is applied to these when --radius changes.
_D = 6.0   # half-separation [a0] in the base configuration (d/R = 0.75)
_B = 2.0

BASE_LOOPS = [
    #   rel. centre [a0]      m direction        n direction
    {"rel": [-_D-_B, -_D+_B, 0.0], "m": [ 1, -1, -1], "n": [ 1,  1,  0]},  # 1
    {"rel": [-_D+_B, +_D+_B, 0.0], "m": [-1, -1,  1], "n": [-1,  1,  0]},  # 2
    {"rel": [+_D-_B, -_D-_B, 0.0], "m": [-1, -1, -1], "n": [-1,  1,  0]},  # 3
    {"rel": [+_D+_B, +_D-_B, 0.0], "m": [ 1, -1,  1], "n": [ 1,  1,  0]},  # 4
]

# Sanity: glide condition b·n = 0 for all base loops
def _dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

for _lp in BASE_LOOPS:
    assert _dot(_lp["m"], _lp["n"]) == 0, \
        f"Glide condition violated: m={_lp['m']}, n={_lp['n']}"

# COM of base quadrupole = average of the 4 centres
_com_check = [sum(lp["rel"][k] for lp in BASE_LOOPS) / 4 for k in range(3)]
assert all(abs(c) < 1e-12 for c in _com_check), \
    f"Base quadrupole COM is not at origin: {_com_check}"


# =============================================================================
# Core builder
# =============================================================================

def build_loops(radius, nrx, nry, nrz, noise_amp, noise_axes, noise_type, rng):
    """
    Return a list of loop dicts and the new box dimensions (in unit cells).

    Each loop dict: {"center": [x,y,z], "m": [...], "n": [...]}
    All positions are in a0 units, origin at box centre.

    Parameters
    ----------
    radius      : float  — new loop radius [a0]
    nrx/y/z     : int    — number of tiles along each axis
    noise_amp   : float  — amplitude of COM noise [a0] (0 = no noise)
    noise_axes  : str    — string of axes to apply noise to, e.g. 'xy', 'xyz'
    noise_type  : str    — 'gaussian' or 'uniform'
    rng         : random.Random instance (or None for no noise)
    """
    scale = radius / BASE_RADIUS   # position scale factor
    cellx  = BASE_CELL_SIZE_X         # a0 — size of one tile along each axis
    celly  = BASE_CELL_SIZE_Y       # a0 — size of one tile along each axis
    cellz  = BASE_CELL_SIZE_Z         # a0 — size of one tile along each axis

    loops = []

    for iz in range(nrz):
        for iy in range(nry):
            for ix in range(nrx):

                # ---- Tile COM in the centred coordinate system ----
                # Centres the whole replicated block at the origin.
                com = [
                    (ix - (nrx - 1) / 2.0) * cellx,
                    (iy - (nry - 1) / 2.0) * celly,
                    (iz - (nrz - 1) / 2.0) * cellz,
                ]

                # ---- Random noise on the tile COM ----
                if noise_amp > 0.0 and rng is not None:
                    for k, axis in enumerate('xyz'):
                        if axis in noise_axes:
                            if noise_type == 'gaussian':
                                com[k] += round(rng.gauss(0.0, noise_amp))
                            else:   # uniform
                                com[k] += round(rng.uniform(-noise_amp, noise_amp))
                                
                # ---- 4 loops per tile ----
                for base in BASE_LOOPS:
                    center = [
                        com[0] + base["rel"][0] * scale,
                        com[1] + base["rel"][1] * scale,
                        com[2] + base["rel"][2] * scale,
                    ]
                    loops.append({
                        "center": center,
                        "m":      base["m"][:],
                        "n":      base["n"][:],
                        "com":    com[:],   # stored for reporting
                    })

    new_box = (BASE_CELL_SIZE_X * nrx,
               BASE_CELL_SIZE_Y * nry,
               BASE_CELL_SIZE_Z * nrz)
    return loops, new_box


# =============================================================================
# Geometry sanity check
# =============================================================================

def check_geometry(loops, box, radius, warn_only=True):
    """Warn (or abort) when a loop extends beyond the periodic box."""
    hx, hy, hz = box[0] / 2.0, box[1] / 2.0, box[2] / 2.0
    bad = False
    for i, lp in enumerate(loops, 1):
        cx, cy, cz = lp["center"]
        checks = [
            (abs(cx) + radius > hx, "x", cx, hx),
            (abs(cy) + radius > hy, "y", cy, hy),
            (abs(cz) + radius > hz, "z", cz, hz),
        ]
        for fails, axis, c, h in checks:
            if fails:
                print(f"  WARNING loop {i}: |{c:.2f}| + R={radius} > {h:.1f} "
                      f"(box half-width along {axis})")
                bad = True
    if bad and not warn_only:
        sys.exit("Geometry error — aborting.")
    return not bad


# =============================================================================
# Input-file writer
# =============================================================================

def write_input(filename, a0, box, radius, loopon_size, loops, fmt):
    """radius and loopon_size in Angstrom; loop centres in a0 (converted here)."""
    nx, ny, nz = box
    with open(filename, "w") as f:
        f.write("general:\n")
        f.write(f"  a0: {a0}\n")
        f.write(f"  Nx: {nx}\n")
        f.write(f"  Ny: {ny}\n")
        f.write(f"  Nz: {nz}\n")
        f.write(f"  Loopon_size: {loopon_size:.10g}\n")
        f.write(f"  Output_format: {fmt}\n")
        f.write("\n")
        f.write("dislocations:\n")
        for i, lp in enumerate(loops, 1):
            cx, cy, cz = lp["center"]
            m = lp["m"]
            n = lp["n"]
            f.write("  - Shape: circular\n")
            f.write(f"    Radius: {radius:.10g}\n")
            f.write(f"    Center: [{cx*a0:.10g}, {cy*a0:.10g}, {cz*a0:.10g}]\n")
            f.write(f"    m: [{m[0]}, {m[1]}, {m[2]}]\n")
            f.write(f"    n: [{n[0]}, {n[1]}, {n[2]}]\n")
    print(f"Written: {filename}")


# =============================================================================
# Main
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="Scale/replicate the quadrupole_circle quadrupole for exec_cuda.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    # -- loop geometry --
    p.add_argument("--radius", type=float, default=BASE_RADIUS * BASE_A0,
                   help=f"Loop radius [Å] (default: {BASE_RADIUS * BASE_A0:g}  — same as base)")
    p.add_argument("--loopon-size", type=float, default=BASE_LOOPON_SIZE * BASE_A0,
                   help=f"Loopon grid step [Å] (default: {BASE_LOOPON_SIZE * BASE_A0:g})")

    # -- replication --
    p.add_argument("--nrx", type=int, default=1,
                   help="Number of quadrupole tiles along x (default: 1)")
    p.add_argument("--nry", type=int, default=1,
                   help="Number of quadrupole tiles along y (default: 1)")
    p.add_argument("--nrz", type=int, default=1,
                   help="Number of quadrupole tiles along z (default: 1)")

    # -- COM noise --
    p.add_argument("--noise-amp", type=float, default=0.0,
                   help="Amplitude of random noise on tile COMs [Å] (default: 0); "
                        "shifts are rounded to whole a0")
    p.add_argument("--noise-type", default="gaussian",
                   choices=["gaussian", "uniform"],
                   help="Noise distribution (default: gaussian)")
    p.add_argument("--noise-axes", default="xyz",
                   help="Axes to apply noise to: any combination of x,y,z "
                        "(default: xyz)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility (default: not set)")

    # -- output --
    p.add_argument("--format", default="lammps",
                   choices=["xyz", "lammps", "both"],
                   help="Output file format (default: lammps)")
    p.add_argument("-o", "--output", default="quadrupole_replicated.yaml",
                   help="Output filename (default: quadrupole_replicated.yaml)")
    p.add_argument("--a0", type=float, default=BASE_A0,
                   help=f"Lattice parameter [Å] (default: {BASE_A0})")

    args = p.parse_args()

    # ---- validate noise axes ----
    noise_axes = args.noise_axes.lower()
    for ch in noise_axes:
        if ch not in 'xyz':
            p.error(f"--noise-axes: invalid character '{ch}', must be a subset of 'xyz'")

    # ---- RNG ----
    rng = None
    if args.noise_amp > 0.0:
        rng = random.Random(args.seed)

    # ---- internal geometry in a0 (command line is in Angstrom) ----
    radius_a0 = args.radius / args.a0
    noise_a0  = args.noise_amp / args.a0

    # ---- build ----
    loops, box = build_loops(
        radius     = radius_a0,
        nrx        = args.nrx,
        nry        = args.nry,
        nrz        = args.nrz,
        noise_amp  = noise_a0,
        noise_axes = noise_axes,
        noise_type = args.noise_type,
        rng        = rng,
    )

    n_tiles = args.nrx * args.nry * args.nrz
    scale   = radius_a0 / BASE_RADIUS
    d_new   = _D * scale   # half-separation after scaling

    # ---- summary ----
    print(f"Base quadrupole  (quadrupole_circle.yaml):")
    print(f"  Radius  = {BASE_RADIUS} a0,  half-sep d = {_D} a0,  "
          f"d/R = {_D/BASE_RADIUS:.3f}")
    print(f"  Box     = {BASE_CELL_SIZE_X} × {BASE_CELL_SIZE_Y} × {BASE_CELL_SIZE_Z} unit cells")
    print()
    print(f"Scaled quadrupole:")
    print(f"  Radius  = {radius_a0:.4g} a0 = {args.radius:.4f} Å   (scale factor = {scale:.4f})")
    print(f"  d       = {d_new:.4f} a0   (half-separation, scaled)")
    print(f"  d/R     = {d_new/radius_a0:.3f}  (preserved)")
    print()
    print(f"Replication:  {args.nrx} × {args.nry} × {args.nrz}  tiles  "
          f"= {n_tiles} quadrupoles  =  {len(loops)} loops")
    print(f"New box:      {box[0]} × {box[1]} × {box[2]} unit cells  "
          f"= {box[0]*args.a0:.2f} × {box[1]*args.a0:.2f} × {box[2]*args.a0:.2f} Å")
    print(f"  ({2*box[0]*box[1]*box[2]} atoms in a perfect BCC cell)")
    print()

    if args.noise_amp > 0.0:
        print(f"COM noise:    {args.noise_type},  amplitude = {args.noise_amp} Å"
              f"  ({noise_a0:.3f} a0)  on axes '{noise_axes}'")
        if args.seed is not None:
            print(f"  seed = {args.seed}")
        print()

    # ---- geometry check ----
    ok = check_geometry(loops, box, radius_a0)
    if ok:
        print(f"Geometry OK:  all loops fit inside the box.")
    print()

    # ---- write ----
    write_input(
        filename   = args.output,
        a0         = args.a0,
        box        = box,
        radius      = args.radius,
        loopon_size = args.loopon_size,
        loops      = loops,
        fmt        = args.format,
    )


if __name__ == "__main__":
    main()
