#!/usr/bin/env python3
"""Write an exec_cuda input for a HEXAGONAL prismatic loop dipole built from six shear loops.

A 1/2<111> Burgers vector b lies in three {110} planes, whose normals are 60 deg
apart around b. Six rectangular loops with the SAME b on +/-n1, +/-n2, +/-n3 are the
faces of a hexagonal tube parallel to b (the glide prism of a BCC prismatic loop).
Every face normal points out of the tube, so the six screw sides (parallel to b)
cancel pairwise between neighbouring faces. What remains are the edge sides at the
two tube ends: two regular hexagonal prismatic loops of opposite sign, in the planes
normal to b, with sides along <112>.

Same construction and checks as make_prismatic_dipole.py (the rhombic, two-plane
version), whose check() and write() are reused: the face planes are put on
mid-planes n.x = +/-(k+1/2) d so exec_cuda's mid-plane snap does not move them, and
Loopon_size is tuned so the effective loop extent (k+1/2)*Loopon_size equals the
hexagon half-side and half the separation. Values actually used are printed.

Usage:
  python3 make_hexagonal_prismatic_dipole.py            # default hexagonal dipole
  python3 make_hexagonal_prismatic_dipole.py --side 25 --separation 100 --b 1 1 1 -o hex2.yaml
"""
import argparse
import itertools
import sys

import numpy as np

import os

from make_prismatic_dipole import OUTDIR, check, mid_spacing, write


def planes_110_containing(b):
    """The three {110} normals n (one sign each) with n.b = 0."""
    found = []
    for n in set(itertools.permutations((1, 1, 0))) | set(itertools.permutations((1, -1, 0))):
        n = np.array(n)
        if n @ b == 0 and not any(np.array_equal(n, -f) or np.array_equal(n, f) for f in found):
            found.append(n)
    return found


def design(a0, b, side, separation, loopon_target):
    b = np.asarray(b, int)
    if sorted(np.abs(b)) != [1, 1, 1]:
        sys.exit(f"b = {b}: give a <111> direction (b = a0/2 [b])")
    normals = planes_110_containing(b)
    assert len(normals) == 3
    bh = b / np.linalg.norm(b)
    d = mid_spacing(normals[0], a0)              # all {110}: same spacing a0/sqrt(2)

    # regular hexagon: apothem s, half-side s / sqrt(3), side 2 s / sqrt(3)
    k_s = max(0, round(side * np.sqrt(3) / 2 / d - 0.5))
    s = (k_s + 0.5) * d                          # face planes on mid-planes
    half_side = s / np.sqrt(3)
    k_p = max(1, round(half_side / loopon_target - 0.5))
    da = half_side / (k_p + 0.5)
    k_m = max(1, round(separation / 2 / da - 0.5))
    half_sep = (k_m + 0.5) * da

    faces = [(sgn * n, s * sgn * n / np.linalg.norm(n))     # centre = apothem point, b.x = 0
             for n in normals for sgn in (1, -1)]
    return dict(script="make_hexagonal_prismatic_dipole.py", b=b, bh=bh, d=d, s=s, k_s=k_s, da=da,
                k_p=k_p, k_m=k_m, half_side=half_side, half_sep=half_sep, faces=faces, normals=normals)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--output", default=os.path.join(OUTDIR, "hexagonal_prismatic_dipole.yaml"))
    p.add_argument("--side", type=float, default=60.67,
                   help="side length of the hexagonal prismatic loops (A); adjusted to mid-planes")
    p.add_argument("--separation", type=float, default=80.3,
                   help="distance between the two prismatic loops along b (A); adjusted to the loopon grid")
    p.add_argument("--loopon-size", type=float, default=1.32,
                   help="target loopon grid step (A); tuned slightly so the edges meet exactly")
    p.add_argument("--cells", type=int, default=50, help="Nx = Ny = Nz (unit cells)")
    p.add_argument("--a0", type=float, default=3.3024)
    p.add_argument("--b", type=int, nargs=3, default=[1, -1, 1], help="Burgers direction <111> (b = a0/2 [b])")
    p.add_argument("--format", default="xyz", choices=["xyz", "lammps", "both", "none"])
    a = p.parse_args()

    D = design(a.a0, a.b, a.side, a.separation, a.loopon_size)
    reach = check(D, a.a0, a.cells)
    write(a.output, D, a.a0, a.cells, a.format)
    print(f"{a.output}: b = a0/2 {[int(x) for x in D['b']]} on {', '.join(str([int(x) for x in n]) for n in D['normals'])}")
    print(f"  hexagon side   {2 * D['half_side']:8.3f} A  (asked {a.side}; face planes at "
          f"n.x = +/-{D['s']:.3f} A = ({D['k_s']}+1/2) d, d = {D['d']:.4f} A)")
    print(f"  separation     {2 * D['half_sep']:8.3f} A  (asked {a.separation})")
    print(f"  Loopon_size    {D['da']:8.6f} A  (asked {a.loopon_size}); "
          f"{(2 * D['k_p'] + 1) * (2 * D['k_m'] + 1)} loopons per face, 6 faces")
    print(f"  box            {a.cells}^3 cells ({2 * a.cells ** 3} atoms); prism reaches {reach:.1f} A "
          f"of the {a.cells * a.a0 / 2:.1f} A half-box")
    print("  checks passed: screw sides cancel pairwise, faces on mid-planes, prism inside the box")


if __name__ == "__main__":
    main()
