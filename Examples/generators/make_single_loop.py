#!/usr/bin/env python3
"""Write exec_cuda inputs for ONE glissile dislocation loop (circular, square, rectangular or
<111> parallelogram).

The slip system is taken from Parameters/BCC.txt (the file DANCE reads), so --system N is
DANCE's own slip-system index: 48 systems, 12 per 1/2<111> Burgers type, each type listing
3 {110}, 3 {112} and 6 {123} planes. b = a0/2 m lies in the loop plane (glissile loop).

Rectangles have their sides along m (screw character, half-length Half_m) and along
p = m x n (edge character, half-length Half_p); --corner-radius rounds the corners.
A parallelogram (MD-like low-temperature BCC shape, {110} planes only) has two screw sides
along b (half-length Half_1 = --half-screw) and two mixed sides along the other <111> of the
plane, at 70.5 deg to b (Half_2 = --half-mixed); exec_cuda uses an oblique loopon grid.
A JSON sidecar (same name, .json) records the ground truth: Burgers type and slip-system
index expected from DANCE, the loop centre after exec_cuda's mid-plane snap, the loopon count.

Usage:
  python3 make_single_loop.py --system 1                                  # circular, R = 8 a0
  python3 make_single_loop.py --shape rectangular --system 16 --corner-radius 13.2
  python3 make_single_loop.py --shape square --family 112 --outdir loops_112/
  python3 make_single_loop.py --family all --units a0 --radius 10         # 48 files
  python3 make_single_loop.py --shape parallelogram --family 110 --outdir para/
  python3 make_single_loop.py --m 1 -1 1 --n 1 1 0 -o my_loop.yaml
"""
import argparse
import itertools
import math
import os
import sys

import exec_input as E

SHAPE_TAG = {"circular": "circ", "square": "sq", "rectangular": "rect", "parallelogram": "para"}


def other_111(m, n):
    """The <111> of plane n other than +/-m (a {110} plane holds exactly two <111>), or None."""
    for v in itertools.product((1, -1), repeat=3):
        if v[0] == 1 and sum(x * y for x, y in zip(v, n)) == 0 and abs(sum(x * y for x, y in zip(v, m))) != 3:
            return list(v)
    return None


def build_loop(a, m, n, unit):
    shape = {"circular": "circular", "parallelogram": "parallelogram"}.get(a.shape, "rectangular")
    loop = dict(shape=shape, center=[c * unit for c in a.center], m=list(m), n=list(n))
    if a.shape == "parallelogram":
        q = other_111(m, n)
        if q is None or sorted(abs(x) for x in m) != [1, 1, 1]:
            sys.exit(f"parallelogram: plane {list(n)} does not contain two <111> with b = {list(m)} "
                     "(only {110} planes do)")
        loop.update(half_1=a.half_screw * unit, half_2=a.half_mixed * unit, dir_1=list(m), dir_2=q)
        return loop
    if a.shape == "circular":
        loop["radius"] = a.radius * unit
    elif a.shape == "square":
        loop.update(half_m=a.half * unit, half_p=a.half * unit, corner_radius=a.corner_radius * unit)
    else:
        loop.update(half_m=a.half_m * unit, half_p=a.half_p * unit, corner_radius=a.corner_radius * unit)
    if loop["shape"] == "rectangular" and loop["corner_radius"] > min(loop["half_m"], loop["half_p"]):
        sys.exit("--corner-radius larger than the smallest half-side (exec_cuda would clamp it)")
    return loop


def expected_content(loop):
    """What DXA + DANCE should find: one closed 1/2<111> loop; nominal screw/edge lengths
    for rectangles (straight parts of the sides, corners excluded)."""
    exp = dict(closed_loops=1, burgers_family="1/2<111>", junctions=0)
    if loop["shape"] == "rectangular":
        rc = loop["corner_radius"]
        exp.update(screw_length_nominal=4 * (loop["half_m"] - rc), edge_length_nominal=4 * (loop["half_p"] - rc))
    if loop["shape"] == "parallelogram":
        exp.update(screw_length_nominal=4 * loop["half_1"], mixed_length_nominal=4 * loop["half_2"],
                   mixed_angle_deg=round(math.degrees(math.acos(1 / 3)), 2))
    return exp


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shape", default="circular", choices=["circular", "square", "rectangular", "parallelogram"])
    sel = p.add_mutually_exclusive_group()
    sel.add_argument("--system", type=int, help="slip-system index 1..48 (BCC.txt order)")
    sel.add_argument("--family", choices=["110", "112", "123", "all"],
                     help="write one file per slip system of this plane family")
    sel.add_argument("--m", type=int, nargs=3, help="Burgers direction, b = a0/2 m (with --n)")
    p.add_argument("--n", type=int, nargs=3, help="loop-plane normal (with --m)")
    p.add_argument("--units", default="A", choices=["A", "a0"], help="unit of every length option")
    p.add_argument("--radius", type=float, help="circular: radius (default 8 a0)")
    p.add_argument("--half", type=float, help="square: half-side (default 8 a0)")
    p.add_argument("--half-m", type=float, help="rectangular: half-length along m, screw sides (default 16 a0)")
    p.add_argument("--half-p", type=float, help="rectangular: half-length along m x n, edge sides (default 8 a0)")
    p.add_argument("--corner-radius", type=float, default=0., help="square/rectangular: 0 = sharp corners")
    p.add_argument("--half-screw", type=float, help="parallelogram: half-length along b, screw sides (default 16 a0)")
    p.add_argument("--half-mixed", type=float,
                   help="parallelogram: half-length along the other <111>, 70.5 deg mixed sides (default 8 a0)")
    p.add_argument("--center", type=float, nargs=3, default=[0., 0., 0.],
                   help="loop centre, origin at box centre (moved onto a mid-plane by exec_cuda)")
    p.add_argument("--loopon-size", type=float, help="loopon grid step (default 0.4 a0 = 1.32096 A)")
    p.add_argument("--cells", type=int, help="Nx = Ny = Nz; default: smallest even cube with a 12 a0 margin "
                   "(the same for every file of one --family run)")
    p.add_argument("--a0", type=float, default=E.A0)
    p.add_argument("--format", default="xyz", choices=["xyz", "lammps", "both", "none"])
    p.add_argument("--bcc", default=E.DEFAULT_BCC, help="BCC.txt used for the slip systems and ground truth")
    p.add_argument("-o", "--output", help="output .yaml (single system only)")
    p.add_argument("--outdir", default=".", help="output folder (default name <shape>_f<family>_sys<NN>.yaml)")
    a = p.parse_args()

    unit = a.a0 if a.units == "a0" else 1.0
    d = lambda default_a0: default_a0 * a.a0 / unit          # defaults are in a0
    a.radius = a.radius if a.radius is not None else d(8.)
    a.half = a.half if a.half is not None else d(8.)
    a.half_m = a.half_m if a.half_m is not None else d(16.)
    a.half_p = a.half_p if a.half_p is not None else d(8.)
    a.half_screw = a.half_screw if a.half_screw is not None else d(16.)
    a.half_mixed = a.half_mixed if a.half_mixed is not None else d(8.)
    da = a.loopon_size * unit if a.loopon_size is not None else 0.4 * a.a0

    bcc = E.load_bcc(a.bcc)
    if a.m is not None:
        if a.n is None:
            sys.exit("--m needs --n")
        if sum(x * y for x, y in zip(a.m, a.n)) != 0:
            sys.exit(f"b = {a.m} is not in plane {a.n} (b.n != 0): not a glissile loop")
        targets = [(a.m, a.n)]
        if E.expected_type_and_sys(bcc, a.m, a.n)[1] is None:
            print(f"WARNING: m {a.m} / n {a.n} is not one of the {len(bcc['systems'])} BCC.txt systems")
    elif a.family is not None:
        targets = [(s["m"], s["n"]) for s in bcc["systems"] if a.family in ("all", s["family"])]
    else:
        idx = a.system or 1
        if not 1 <= idx <= len(bcc["systems"]):
            sys.exit(f"--system must be in 1..{len(bcc['systems'])}")
        targets = [(bcc["systems"][idx - 1]["m"], bcc["systems"][idx - 1]["n"])]
    if a.output and len(targets) > 1:
        sys.exit("-o is for a single system; use --outdir with --family")

    loops = [build_loop(a, m, n, unit) for m, n in targets]
    cells = a.cells or max(E.auto_cells([l], a.a0) for l in loops)
    general = E.general_block(a.a0, cells, da, a.format)
    os.makedirs(a.outdir, exist_ok=True)

    for loop in loops:
        E.check_in_box([loop], a.a0, (cells,) * 3)
        ty, idx = E.expected_type_and_sys(bcc, loop["m"], loop["n"])
        fam = E.plane_family(loop["n"])
        name = a.output or os.path.join(
            a.outdir, f"{SHAPE_TAG[a.shape]}{'_rounded' if a.corner_radius > 0 else ''}_f{fam}_"
                      + (f"sys{idx:02d}" if idx else f"m{''.join(map(str, loop['m']))}_n{''.join(map(str, loop['n']))}")
                      + ".yaml")
        m_s, n_s = " ".join(map(str, loop["m"])), " ".join(map(str, loop["n"]))
        header = [f"Single glissile {a.shape} loop, generated by make_single_loop.py",
                  f"b = a0/2 [{m_s}] on ({n_s}) {{{fam}}}: DANCE Burgers type {ty}, slip system {idx}"]
        E.write_yaml(name, general, [loop], header)
        E.write_sidecar(os.path.splitext(name)[0] + ".json", bcc, general, [loop], expected_content(loop))
        print(f"{name}: sys {idx} type {ty} {{{fam}}}  m [{m_s}] n [{n_s}]  "
              f"{E.count_loopons(loop, da)} loopons, box {cells}^3")


if __name__ == "__main__":
    main()
