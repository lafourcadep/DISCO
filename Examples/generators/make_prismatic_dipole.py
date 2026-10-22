#!/usr/bin/env python3
"""Write an exec_cuda input for a prismatic loop dipole built from four shear loops.

Four rectangular loops carry the SAME Burgers vector b and lie on the two slip
planes n1, n2 that contain b. They are the faces of a tube parallel to b (rhombic
cross-section). Every face normal points out of the tube, so the tube is one
consistently oriented surface: the screw sides (parallel to b) of neighbouring
faces coincide with opposite line senses and cancel. What remains are the edge
sides at the two tube ends, normal to b: two rhombic prismatic loops of opposite
sign (a prismatic dipole), separated by --separation along b.

Exact cancellation needs the shared edges to coincide, so the requested sizes are
adjusted and the values actually used are printed:
  * the face planes are put on mid-planes between atomic planes, n.x = +/-(k+1/2) d,
    so exec_cuda's automatic mid-plane snap does not move them (sets the side length);
  * Loopon_size is tuned (close to --loopon-size) so that the effective loop extent,
    (k+1/2) * Loopon_size, equals the rhombus half-side and half the separation.

Usage:
  python3 make_prismatic_dipole.py                       # writes ./prismatic_dipole.yaml (+ .json)
  python3 make_prismatic_dipole.py --side 40 --separation 120 --cells 60 -o dipole2.yaml
  python3 make_prismatic_dipole.py --b 1 1 1 --n1 1 -1 0 --n2 0 1 -1
"""
import argparse
import os
import sys

import numpy as np

import exec_input as E

OUTDIR = "."                        # default output folder


def mid_spacing(n, a0):
    """Spacing of the atomic planes of family n in BCC (as snapToMidPlane in loop.cpp)."""
    n = np.asarray(n)
    return a0 / np.linalg.norm(n) if n.sum() % 2 == 0 else a0 / (2 * np.linalg.norm(n))


def design(a0, b, n1, n2, side, separation, loopon_target):
    b, n1, n2 = (np.asarray(v, int) for v in (b, n1, n2))
    for n in (n1, n2):
        if b @ n != 0:
            sys.exit(f"b = {b} does not lie in plane {n} (b.n = {b @ n})")
    if abs(np.linalg.norm(np.cross(n1, n2))) < 1e-12:
        sys.exit("n1 and n2 must be different planes")
    d = mid_spacing(n1, a0)
    if abs(d - mid_spacing(n2, a0)) > 1e-12:
        sys.exit("n1 and n2 must have the same plane spacing (e.g. both {110})")

    bh = b / np.linalg.norm(b)
    nh1, nh2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    sin = abs(nh2 @ np.cross(bh, nh1))          # rhombus: half-side = s / sin
    k_s = max(0, round(side / 2 * sin / d - 0.5))
    s = (k_s + 0.5) * d                          # face planes on mid-planes
    half_side = s / sin
    k_p = max(1, round(half_side / loopon_target - 0.5))
    da = half_side / (k_p + 0.5)                 # effective extent (k+1/2) da = half-side
    k_m = max(1, round(separation / 2 / da - 0.5))
    half_sep = (k_m + 0.5) * da

    faces = []                                   # (integer outward normal, centre)
    for n_int, other in ((n1, nh2), (n2, nh1)):
        for sgn in (1, -1):
            nh = sgn * n_int / np.linalg.norm(n_int)
            c = np.linalg.solve([nh, other, bh], [s, 0., 0.])   # mid-side, b.x = 0
            faces.append((sgn * n_int, c))
    return dict(script="make_prismatic_dipole.py", b=b, bh=bh, d=d, s=s, k_s=k_s, da=da, k_p=k_p, k_m=k_m,
                half_side=half_side, half_sep=half_sep, faces=faces)


def check(D, a0, cells):
    """Screw sides cancel pairwise, ends are normal to b, no snap shift, prism in the box."""
    bh, E_m, E_p = D["bh"], D["half_sep"], D["half_side"]
    segs = []
    for n_int, c in D["faces"]:
        nh = n_int / np.linalg.norm(n_int)
        p = np.cross(bh, nh)
        corners = [c + u * E_m * bh + v * E_p * p for u, v in ((1, 1), (-1, 1), (-1, -1), (1, -1))]
        if np.cross(corners[1] - corners[0], corners[2] - corners[1]) @ nh < 0:
            corners = corners[::-1]              # right-hand line sense about n
        segs += [(corners[i], corners[(i + 1) % 4]) for i in range(4)]
    screw = [sg for sg in segs if abs(abs((sg[1] - sg[0]) @ bh) - 2 * E_m) < 1e-6]
    edge = [sg for sg in segs if abs((sg[1] - sg[0]) @ bh) < 1e-6]
    assert len(screw) == len(edge) == 2 * len(D["faces"])
    assert all(any(np.allclose(a, b2) and np.allclose(b, a2) for a2, b2 in screw) for a, b in screw), \
        "screw sides do not cancel"
    for n_int, c in D["faces"]:                  # same formula as snapToMidPlane in loop.cpp
        nh = n_int / np.linalg.norm(n_int)
        sp = mid_spacing(n_int, a0)
        s0 = -0.5 * a0 * cells * nh.sum()
        sc = c @ nh
        assert abs(s0 + (np.floor((sc - s0) / sp + 1e-9) + 0.5) * sp - sc) < 1e-9, "face not on a mid-plane"
    reach = max(np.abs(p).max() for sg in segs for p in sg)
    half_box = cells * a0 / 2
    if reach > half_box - 10.:
        sys.exit(f"prism reaches {reach:.1f} A but the half-box is {half_box:.1f} A: "
                 f"increase --cells (need > {2 * (reach + 10.) / a0:.0f})")
    return reach


def write(path, D, a0, cells, fmt):
    g = lambda x: f"{x:.10g}"
    b = D["b"]
    lines = [f"# Prismatic loop dipole from {len(D['faces'])} shear loops; generated by {D['script']}",
             f"# b = a0/2 [{' '.join(map(str, b))}], loop side {2 * D['half_side']:.2f} A, "
             f"loops at b.x = +/-{D['half_sep']:.2f} A",
             "general:", f"  a0: {a0}", f"  Nx: {cells}", f"  Ny: {cells}", f"  Nz: {cells}",
             f"  Loopon_size: {g(D['da'])}", f"  Output_format: {fmt}", "", "dislocations:"]
    for n_int, c in D["faces"]:
        lines += ["  - Shape: rectangular",
                  f"    Half_m: {g((D['k_m'] + 0.25) * D['da'])}",   # quarter step past the last
                  f"    Half_p: {g((D['k_p'] + 0.25) * D['da'])}",   # grid row: no rounding loss
                  "    Corner_radius: 0",
                  f"    Center: [{', '.join(g(x) for x in c)}]",
                  f"    m: [{', '.join(map(str, b))}]",
                  f"    n: [{', '.join(str(int(x)) for x in n_int)}]"]
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    open(path, "w").write("\n".join(lines) + "\n")
    write_sidecar(path, D, a0, cells, fmt)


def write_sidecar(path, D, a0, cells, fmt):
    """JSON ground truth next to the YAML (exec_input.write_sidecar): after the screw sides
    cancel, two closed prismatic 1/2<111> loops of opposite sign, normal to b, 0 screw."""
    loops = [dict(shape="rectangular", half_m=(D["k_m"] + 0.25) * D["da"], half_p=(D["k_p"] + 0.25) * D["da"],
                  corner_radius=0., center=list(c), m=[int(x) for x in D["b"]], n=[int(x) for x in n_int])
             for n_int, c in D["faces"]]
    general = E.general_block(a0, cells, D["da"], fmt)
    nsides = len(D["faces"])
    expected = dict(closed_loops=2, burgers_family="1/2<111>", junctions=0, screw_length=0.,
                    prismatic_loops=dict(shape="rhombus" if nsides == 4 else "hexagon", sides=nsides,
                                         side_length=2 * D["half_side"], normal=[int(x) for x in D["b"]],
                                         positions_along_b=[-D["half_sep"], D["half_sep"]],
                                         signs="opposite (prismatic dipole)"))
    E.write_sidecar(os.path.splitext(path)[0] + ".json", E.load_bcc(), general, loops, expected,
                    extra=dict(note="shear loops of the same b on the faces of a glide prism; the screw sides "
                                    "cancel pairwise, leaving the two end loops"), script=D["script"])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--output", default=os.path.join(OUTDIR, "prismatic_dipole.yaml"))
    p.add_argument("--side", type=float, default=56.6,
                   help="side length of the rhombic prismatic loops (A); adjusted to mid-planes")
    p.add_argument("--separation", type=float, default=80.3,
                   help="distance between the two prismatic loops along b (A); adjusted to the loopon grid")
    p.add_argument("--loopon-size", type=float, default=1.32,
                   help="target loopon grid step (A); tuned slightly so the edges meet exactly")
    p.add_argument("--cells", type=int, default=50, help="Nx = Ny = Nz (unit cells)")
    p.add_argument("--a0", type=float, default=3.3024)
    p.add_argument("--b", type=int, nargs=3, default=[1, -1, 1], help="Burgers direction (b = a0/2 [b])")
    p.add_argument("--n1", type=int, nargs=3, default=[1, 1, 0], help="first slip plane containing b")
    p.add_argument("--n2", type=int, nargs=3, default=[0, 1, 1], help="second slip plane containing b")
    p.add_argument("--format", default="xyz", choices=["xyz", "lammps", "both", "none"])
    a = p.parse_args()

    D = design(a.a0, a.b, a.n1, a.n2, a.side, a.separation, a.loopon_size)
    reach = check(D, a.a0, a.cells)
    write(a.output, D, a.a0, a.cells, a.format)
    print(f"{a.output}: b = a0/2 {[int(x) for x in D['b']]} on {list(a.n1)} and {list(a.n2)}")
    print(f"  loop side      {2 * D['half_side']:8.3f} A  (asked {a.side}; face planes at "
          f"n.x = +/-{D['s']:.3f} A = ({D['k_s']}+1/2) d, d = {D['d']:.4f} A)")
    print(f"  separation     {2 * D['half_sep']:8.3f} A  (asked {a.separation})")
    print(f"  Loopon_size    {D['da']:8.6f} A  (asked {a.loopon_size}); "
          f"{(2 * D['k_p'] + 1) * (2 * D['k_m'] + 1)} loopons per face")
    print(f"  box            {a.cells}^3 cells ({2 * a.cells ** 3} atoms); prism reaches {reach:.1f} A "
          f"of the {a.cells * a.a0 / 2:.1f} A half-box")
    print("  checks passed: screw sides cancel pairwise, faces on mid-planes, prism inside the box")


if __name__ == "__main__":
    main()
