#!/usr/bin/env python3
"""
generate_quadrupole.py
======================
Generate an input file (quadrupole of circular loops) for the Kroupa CUDA code.

Two quadrupole modes
--------------------
junction (default)
    Two BCC {110}<111> slip systems whose Burgers vectors satisfy Frank's rule:
        a0/2 [1,-1, 1](110)  +  a0/2 [1, 1,-1](011)  →  a0 [1,0,0]
        |b_junction|² = a0²  <  |b1|² + |b2|² = 3a0²/2   ✓
    The junction product a0[1,0,0] is a hard (Hirth-lock type) dislocation:
    it cannot dissociate back into a0/2<111> partials — effective obstacle.
    Layout (centres in a0, origin at box centre):
        Loop 1:  +b_A  at (−d, −d, 0)   system A (110)
        Loop 2:  −b_A  at (+d, +d, 0)   system A (110)  — stress partner
        Loop 3:  +b_B  at (+d, −d, 0)   system B (011)
        Loop 4:  −b_B  at (−d, +d, 0)   system B (011)  — stress partner
    Loops 1 & 3 (and 2 & 4) are the closest cross-system pairs (distance 2d);
    they interact along the intersection line of their glide planes and can
    form a sessile junction.

annihilation
    Two dipoles, each on a different {110} glide plane, with opposite Burgers
    vectors.  Loops within each dipole are co-planar (same glide plane) and
    carry opposite b; they annihilate when they meet upon expansion.
        Plane  (110) [x+y = 0]:  b = a0/2 [−1, 1, 1]
            Loop 1:  +b  at (−d, +d, 0)
            Loop 2:  −b  at (+d, −d, 0)
        Plane (1-10) [x−y = 0]:  b = a0/2 [ 1, 1, 1]
            Loop 3:  +b  at (−d, −d, 0)
            Loop 4:  −b  at (+d, +d, 0)
    Note: a cross-system junction reaction is also energetically possible
    (b1 − b3 = a0[−1,0,0]), but the geometry favours intra-dipole annihilation.

Usage
-----
    python generate_quadrupole.py [options]

    python generate_quadrupole.py --mode junction  --radius 33.024 -o quad_junc.yaml
    python generate_quadrupole.py --mode annihilation --radius 26.4192 -o quad_anni.yaml

All lengths on the command line (--radius, --distance, --loopon-size) are in
Angstrom, like in the exec_cuda input file. Loop layouts below are in a0.
    python generate_quadrupole.py --list-slip-systems
"""

import argparse
import math
import sys

# =============================================================================
# BCC {110}<111> slip systems — complete 12-system table
# Convention: b and n are integer Miller indices (not normalised).
# Glide condition b·n = 0 is verified below.
# =============================================================================
BCC_SLIP_SYSTEMS = [
    # plane (110)  n = [1, 1, 0]
    ([ 1, -1,  1], [ 1,  1,  0]),   #  1
    ([-1,  1,  1], [ 1,  1,  0]),   #  2
    # plane (1-10) n = [1,-1, 0]
    ([ 1,  1,  1], [ 1, -1,  0]),   #  3
    ([-1, -1,  1], [ 1, -1,  0]),   #  4
    # plane (101)  n = [1, 0, 1]
    ([-1,  1,  1], [ 1,  0,  1]),   #  5
    ([ 1,  1, -1], [ 1,  0,  1]),   #  6
    # plane (10-1) n = [1, 0,-1]
    ([ 1,  1,  1], [ 1,  0, -1]),   #  7
    ([ 1, -1,  1], [ 1,  0, -1]),   #  8
    # plane (011)  n = [0, 1, 1]
    ([ 1, -1,  1], [ 0,  1,  1]),   #  9
    ([ 1,  1, -1], [ 0,  1,  1]),   # 10
    # plane (01-1) n = [0, 1,-1]
    ([ 1,  1,  1], [ 0,  1, -1]),   # 11
    ([-1,  1,  1], [ 0,  1, -1]),   # 12
]

# Quick sanity check (runs once at import)
for _b, _n in BCC_SLIP_SYSTEMS:
    assert _b[0]*_n[0] + _b[1]*_n[1] + _b[2]*_n[2] == 0, \
        f"Glide condition violated: b={_b}, n={_n}"


# =============================================================================
# Helpers
# =============================================================================

def neg(v):
    return [-v[0], -v[1], -v[2]]

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def frank(b1, b2):
    """Return (|b1+b2|², |b1|²+|b2|²) in units of (a0/2)²."""
    b3 = [b1[i] + b2[i] for i in range(3)]
    return dot(b3, b3), dot(b1, b1) + dot(b2, b2)


# =============================================================================
# Quadrupole configurations
# =============================================================================

def quadrupole_junction(d):
    """
    System A: b = a0/2 [1,-1, 1]  n = (110)   b·n = 0 ✓
    System B: b = a0/2 [1, 1,-1]  n = (011)   b·n = 0 ✓
    Junction: b_A + b_B = a0/2 [2,0,0] = a0 [1,0,0]
    Frank:    |[2,0,0]|² = 4  <  |[1,-1,1]|² + |[1,1,-1]|² = 6   ✓
    """
    bA, nA = [ 1, -1,  1], [1, 1, 0]
    bB, nB = [ 1,  1, -1], [0, 1, 1]
    return [
        {"b": bA,      "n": nA, "center": [-d, -d, 0]},  # +b_A
        {"b": neg(bA), "n": nA, "center": [ d,  d, 0]},  # −b_A  (stress partner)
        {"b": bB,      "n": nB, "center": [ d, -d, 0]},  # +b_B
        {"b": neg(bB), "n": nB, "center": [-d,  d, 0]},  # −b_B  (stress partner)
    ]


def quadrupole_annihilation(d):
    """
    Plane (110)  b = a0/2 [-1, 1, 1]:
        Loop 1 centre (−d, +d, 0)  →  x+y = 0   ✓
        Loop 2 centre (+d, −d, 0)  →  x+y = 0   ✓  (opposite b → annihilate)
    Plane (1-10) b = a0/2 [1, 1, 1]:
        Loop 3 centre (−d, −d, 0)  →  x−y = 0   ✓
        Loop 4 centre (+d, +d, 0)  →  x−y = 0   ✓  (opposite b → annihilate)
    """
    bA, nA = [-1,  1,  1], [ 1,  1,  0]   # plane (110),   x+y = 0
    bB, nB = [ 1,  1,  1], [ 1, -1,  0]   # plane (1-10),  x-y = 0
    return [
        {"b": bA,      "n": nA, "center": [-d,  d, 0]},  # plane 1, +b
        {"b": neg(bA), "n": nA, "center": [ d, -d, 0]},  # plane 1, −b
        {"b": bB,      "n": nB, "center": [-d, -d, 0]},  # plane 2, +b
        {"b": neg(bB), "n": nB, "center": [ d,  d, 0]},  # plane 2, −b
    ]


# =============================================================================
# Input-file writer
# =============================================================================

def write_input(filename, a0, nx, ny, nz, radius, loopon_size, loops, fmt):
    """radius and loopon_size in Angstrom; loop centres in a0 (converted here)."""
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
            b, n, c = lp["b"], lp["n"], lp["center"]
            f.write("  - Shape: circular\n")
            f.write(f"    Radius: {radius:.10g}\n")
            f.write(f"    Center: [{c[0]*a0:.10g}, {c[1]*a0:.10g}, {c[2]*a0:.10g}]\n")
            f.write(f"    m: [{b[0]}, {b[1]}, {b[2]}]\n")
            f.write(f"    n: [{n[0]}, {n[1]}, {n[2]}]\n")
    print(f"Written: {filename}")


# =============================================================================
# Main
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="Generate a quadrupole loop input for the Kroupa CUDA code.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    p.add_argument("-o", "--output",      default="quadrupole.yaml",
                   help="Output filename (default: quadrupole.yaml)")
    p.add_argument("--mode",              default="junction",
                   choices=["junction", "annihilation"],
                   help="Quadrupole type (default: junction)")
    p.add_argument("--a0",                type=float, default=3.3024,
                   help="Lattice parameter [Å] (default: 3.3024)")
    p.add_argument("--nx",                type=int,   default=40,
                   help="Box size along x [unit cells] (default: 40)")
    p.add_argument("--ny",                type=int,   default=40,
                   help="Box size along y [unit cells] (default: 40)")
    p.add_argument("--nz",                type=int,   default=40,
                   help="Box size along z [unit cells] (default: 40)")
    p.add_argument("--radius",            type=float, default=33.024,
                   help="Loop radius [Å] (default: 33.024 = 10 a0)")
    p.add_argument("--loopon-size",       type=float, default=1.32096,
                   help="Loopon grid step [Å] (default: 1.32096 = 0.4 a0)")
    p.add_argument("--distance",          type=float, default=None,
                   help="Half-separation of loop centres [Å] "
                        "(default: radius/2)")
    p.add_argument("--format",            default="both",
                   choices=["xyz", "lammps", "both"],
                   help="Output file format (default: both)")
    p.add_argument("--list-slip-systems", action="store_true",
                   help="Print all 12 BCC {110}<111> slip systems and exit")

    args = p.parse_args()

    # ---- optional: print slip-system table ----
    if args.list_slip_systems:
        print("BCC {110}<111> slip systems  (b·n = 0 for all):\n")
        print(f"{'#':>3}  {'b  [Miller × a0/2]':>20}  {'n  [Miller]':>16}  b·n")
        print("-" * 50)
        for i, (b, n) in enumerate(BCC_SLIP_SYSTEMS, 1):
            print(f"{i:>3}  {str(b):>20}  {str(n):>16}  {dot(b,n)}")
        sys.exit(0)

    # ---- internal geometry in a0 (command line is in Angstrom) ----
    radius_a0 = args.radius / args.a0
    # default distance: radius / 2  (loops overlap slightly)
    d = (args.distance if args.distance is not None else args.radius * 0.5) / args.a0

    # ---- build loop descriptors ----
    if args.mode == "junction":
        loops = quadrupole_junction(d)
        bA, bB = [1,-1,1], [1,1,-1]
        ep, es = frank(bA, bB)
        b3 = [bA[i]+bB[i] for i in range(3)]

        print("Mode: junction")
        print(f"  System A:  b = a0/2 {bA}   n = (110)")
        print(f"  System B:  b = a0/2 {bB}   n = (011)")
        print(f"  Junction reaction:  a0/2{bA} + a0/2{bB} → a0/2{b3} = a0[1,0,0]")
        print(f"  Frank's rule:  |b_A+b_B|² = {ep} (a0/2)²  <  "
              f"|b_A|²+|b_B|² = {es} (a0/2)²  ✓  (junction favourable)")
        print(f"  a0[1,0,0] is a hard (Hirth-lock type) junction: "
              f"energetically stable obstacle.")
        print()
        print(f"  Loop centres [a0 from box centre]:")
        for i, lp in enumerate(loops, 1):
            c, b = lp["center"], lp["b"]
            sys_label = "A" if i <= 2 else "B"
            sign = "+" if b[0] > 0 else "−"
            print(f"    Loop {i}  {sign}b_{sys_label}  n={lp['n']}  "
                  f"centre = ({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f})")
        print()
        print(f"  Cross-system pairs that can form junctions:")
        print(f"    Loops 1 & 3  (+b_A and +b_B)  distance = {2*d:.2f} a0")
        print(f"    Loops 2 & 4  (−b_A and −b_B)  distance = {2*d:.2f} a0")

    else:  # annihilation
        loops = quadrupole_annihilation(d)
        bA, bB = [-1,1,1], [1,1,1]

        # verify co-planarity for each dipole
        assert loops[0]["center"][0] + loops[0]["center"][1] == \
               loops[1]["center"][0] + loops[1]["center"][1], \
            "Annihilation dipole 1: loops not co-planar!"
        assert loops[2]["center"][0] - loops[2]["center"][1] == \
               loops[3]["center"][0] - loops[3]["center"][1], \
            "Annihilation dipole 2: loops not co-planar!"

        print("Mode: annihilation")
        print(f"  Dipole 1 — plane (110)  b = a0/2 {bA}")
        print(f"    Loop 1: +b  centre ({-d:.1f}, {d:.1f}, 0)   glide plane x+y = 0")
        print(f"    Loop 2: −b  centre ({d:.1f}, {-d:.1f}, 0)   glide plane x+y = 0")
        print(f"  Dipole 2 — plane (1-10) b = a0/2 {bB}")
        print(f"    Loop 3: +b  centre ({-d:.1f}, {-d:.1f}, 0)  glide plane x−y = 0")
        print(f"    Loop 4: −b  centre ({d:.1f}, {d:.1f}, 0)   glide plane x−y = 0")
        print()
        print(f"  Loops 1+2 are co-planar with opposite b  →  annihilate when they meet.")
        print(f"  Loops 3+4 are co-planar with opposite b  →  annihilate when they meet.")

        # note on cross-system reaction
        ep, es = frank(bA, neg(bB))
        if ep < es:
            b3 = [bA[i] - bB[i] for i in range(3)]
            print()
            print(f"  Note: cross-system junction {bA}+(−{bB})={b3} also satisfies Frank's rule")
            print(f"    ({ep}(a0/2)² < {es}(a0/2)²) but intra-dipole annihilation dominates")
            print(f"    for the chosen geometry.")

    print()
    print(f"Box:         {args.nx} × {args.ny} × {args.nz}  unit cells")
    print(f"a0:          {args.a0} Å")
    print(f"Loop radius: {radius_a0:.4g} a0  =  {args.radius:.2f} Å")
    print(f"Separation:  d = {d} a0  =  {d * args.a0:.2f} Å  [centre to origin]")
    print()

    # ---- geometry sanity check ----
    half = [args.nx / 2.0, args.ny / 2.0, args.nz / 2.0]
    warned = False
    for i, lp in enumerate(loops, 1):
        for k, (ck, hk) in enumerate(zip(lp["center"], half)):
            if abs(ck) + radius_a0 > hk:
                axis = "xyz"[k]
                print(f"  WARNING loop {i}: extends outside box along {axis}  "
                      f"(|{ck}| + R={radius_a0:.4g} a0 > N{axis}/2={hk})")
                warned = True
    if warned:
        print()

    write_input(
        args.output,
        a0=args.a0, nx=args.nx, ny=args.ny, nz=args.nz,
        radius=args.radius, loopon_size=args.loopon_size,
        loops=loops, fmt=args.format)


if __name__ == "__main__":
    main()
