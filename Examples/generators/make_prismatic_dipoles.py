#!/usr/bin/env python3
"""Write an exec_cuda input with N prismatic loop dipoles placed at random in the box.

Each dipole is built as in make_prismatic_dipole.py (rhombic: 4 shear loops on the two
{110} planes of one pair containing b) or make_hexagonal_prismatic_dipole.py (hexagonal:
6 shear loops on the three {110} planes containing b): the screw sides cancel pairwise and
leave two prismatic loops of opposite sign, normal to b, separated along b.

Per dipole, at random (--seed):
  b             one of the eight 1/2<111> vectors (or the --b direction for all; its sign
                still flips with --random-sign), so the loop of positive sign is at +b or -b;
  planes        rhombic: one of the three {110} pairs containing b;
  separation    uniform in [--separation MIN MAX], on the loopon grid ((k + 1/2) Loopon_size);
  position      a random lattice vector (whole a0 on each axis), so every face plane
                stays on its mid-plane and exec_cuda's snap does not move it.
A dipole is rejected and redrawn when its bounding sphere comes closer than --gap to that of
another dipole (periodic minimum image).

One sample holds one type: the loopon step that makes the screw sides cancel exactly is
set by the side length, and no step fits a rhombic and a hexagonal dipole at once (their
half-sides are odd multiples of d/sqrt(3) and of d/(2 sqrt(3)), d = {110} spacing). All
dipoles of a sample therefore have the same side; the separation can differ.

Usage:
  python3 make_prismatic_dipoles.py --count 10                       # 10 hexagonal dipoles, random b
  python3 make_prismatic_dipoles.py --count 20 --type rhombic --side 40 --separation 60 120 --seed 3
  python3 make_prismatic_dipoles.py --count 8 --b 1 1 1 --cells 120 --format lammps
"""
import argparse
import itertools
import math
import os
import sys

import numpy as np

import exec_input as E
import make_hexagonal_prismatic_dipole as HEX
import make_prismatic_dipole as RHO

B111 = [np.array(v) for v in itertools.product((1, -1), repeat=3)]


def design_one(a, b, pair, separation):
    if a.type == "hexagonal":
        return HEX.design(a.a0, b, a.side, separation, a.loopon_size)
    return RHO.design(a.a0, b, pair[0], pair[1], a.side, separation, a.loopon_size)


def face_loops(D):
    """The shear loops of one centred dipole, as exec_input loop dicts."""
    return [dict(shape="rectangular", half_m=(D["k_m"] + 0.25) * D["da"], half_p=(D["k_p"] + 0.25) * D["da"],
                 corner_radius=0., center=np.asarray(c, float), m=[int(x) for x in D["b"]],
                 n=[int(x) for x in n_int]) for n_int, c in D["faces"]]


def bounding_radius(D):
    """Largest distance from the dipole centre to a corner of its shear loops."""
    r = 0.
    for n_int, c in D["faces"]:
        p = np.cross(D["bh"], n_int / np.linalg.norm(n_int))
        for u, v in itertools.product((1, -1), repeat=2):
            r = max(r, np.linalg.norm(c + u * D["half_sep"] * D["bh"] + v * D["half_side"] * p))
    return r


def min_image(d, box):
    return d - box * np.round(d / box)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--count", type=int, required=True, help="number of dipoles")
    p.add_argument("--type", default="hexagonal", choices=["hexagonal", "rhombic"])
    p.add_argument("--side", type=float, help="side length of the prismatic loops (A; default 60.67 hexagonal, "
                   "56.6 rhombic); adjusted to mid-planes, the same for every dipole")
    p.add_argument("--separation", type=float, nargs="+", default=[80.3], metavar="A",
                   help="distance between the two prismatic loops along b (A): one value, or MIN MAX for a "
                        "uniform random value per dipole; adjusted to the loopon grid")
    p.add_argument("--b", type=int, nargs=3, help="Burgers direction <111> for every dipole "
                   "(default: random among the four <111>)")
    p.add_argument("--random-sign", action=argparse.BooleanOptionalAction, default=True,
                   help="flip the sign of b at random (which end loop is positive); default on")
    p.add_argument("--cells", type=int, help="Nx = Ny = Nz; default: a cube where the dipoles' bounding "
                   "spheres fill --fill of the volume")
    p.add_argument("--fill", type=float, default=0.15,
                   help="default box: volume fraction of the bounding spheres (with --gap)")
    p.add_argument("--gap", type=float, default=10.0, help="minimum distance between the bounding spheres "
                   "of two dipoles (A, periodic minimum image)")
    p.add_argument("--loopon-size", type=float, default=1.32,
                   help="target loopon grid step (A); tuned slightly so the edges meet exactly")
    p.add_argument("--a0", type=float, default=E.A0)
    p.add_argument("--format", default="xyz", choices=["xyz", "lammps", "both", "none"])
    p.add_argument("--seed", type=int, help="random seed")
    p.add_argument("--tries", type=int, default=10000, help="placement attempts per dipole")
    p.add_argument("-o", "--output", help="output .yaml (default: prismatic_dipoles_<count>_<type>.yaml)")
    p.add_argument("--outdir", default=".")
    a = p.parse_args()

    if a.count < 1:
        sys.exit("--count must be >= 1")
    if len(a.separation) not in (1, 2) or min(a.separation) <= 0:
        sys.exit("--separation: one positive value, or MIN MAX")
    if len(a.separation) == 2 and a.separation[0] > a.separation[1]:
        sys.exit("--separation MIN MAX: MIN > MAX")
    if a.b is not None and sorted(abs(x) for x in a.b) != [1, 1, 1]:
        sys.exit(f"--b {a.b}: give a <111> direction")
    if a.side is None:
        a.side = 60.67 if a.type == "hexagonal" else 56.6
    rng = np.random.default_rng(a.seed)

    # ---- draw the dipoles (b, planes, separation), centred on the origin ----
    dirs = [np.array(a.b)] if a.b is not None else [v for v in B111 if v[0] > 0]   # four <111>
    dipoles = []
    for _ in range(a.count):
        b = dirs[int(rng.integers(len(dirs)))]
        if a.random_sign and rng.random() < 0.5:
            b = -b
        pair = None
        if a.type == "rhombic":
            normals = HEX.planes_110_containing(b)
            i, j = sorted(rng.choice(3, size=2, replace=False))
            pair = (normals[i], normals[j])
        sep = a.separation[0] if len(a.separation) == 1 else float(rng.uniform(*a.separation))
        D = design_one(a, b, pair, sep)
        RHO.check(D, a.a0, 2 * math.ceil(bounding_radius(D) / a.a0) + 10)   # screw sides cancel, mid-planes
        dipoles.append(D)
    das = sorted({round(float(D["da"]), 10) for D in dipoles})
    if len(das) > 1:
        sys.exit(f"internal error: dipoles with different Loopon_size {das}")
    da = dipoles[0]["da"]
    radii = [bounding_radius(D) for D in dipoles]

    # ---- box ----
    if a.cells:
        cells = a.cells
    else:
        vol = sum(4 / 3 * math.pi * (r + a.gap / 2) ** 3 for r in radii) / a.fill
        cells = math.ceil(vol ** (1 / 3) / a.a0)
        cells = max(cells, math.ceil(2 * (max(radii) + a.gap) / a.a0))
    box = cells * a.a0
    if 2 * (max(radii) + a.gap / 2) > box:
        sys.exit(f"a dipole spans {2 * max(radii):.1f} A, more than the {box:.1f} A box: increase --cells")

    # ---- random placement on lattice vectors, no overlap of bounding spheres ----
    order = np.argsort(radii)[::-1]                           # large ones first
    centres = [None] * a.count
    for t in order:
        for _ in range(a.tries):
            c = rng.integers(-(cells // 2), cells - cells // 2, size=3) * a.a0
            if all(np.linalg.norm(min_image(c - centres[u], box)) >= radii[t] + radii[u] + a.gap
                   for u in range(a.count) if centres[u] is not None):
                centres[t] = c.astype(float)
                break
        else:
            sys.exit(f"could not place dipole {t + 1} after {a.tries} tries: increase --cells, "
                     f"lower --fill or --gap, or reduce --count")

    # ---- loops, checks, records ----
    loops, records = [], []
    for t, (D, c) in enumerate(zip(dipoles, centres)):
        first = len(loops) + 1
        for l in face_loops(D):
            l["center"] = l["center"] + c
            moved = np.linalg.norm(E.snap_center(l["center"], l["n"], a.a0, (cells,) * 3) - l["center"])
            if moved > 1e-6:
                sys.exit(f"internal error: dipole {t + 1}: a face is not on a mid-plane (moved {moved:.2e} A)")
            loops.append(l)
        bh = D["bh"]
        records.append(dict(dipole=t + 1, shear_loops=list(range(first, len(loops) + 1)),
                            b=[int(x) for x in D["b"]], centre=c,
                            planes=[[int(x) for x in n] for n, _ in D["faces"][::2]],
                            side_length=2 * D["half_side"], separation=2 * D["half_sep"],
                            prismatic_loop_centres=[c + D["half_sep"] * bh, c - D["half_sep"] * bh],
                            bounding_radius=radii[t]))
    general = E.general_block(a.a0, cells, da, a.format)
    nsides = 6 if a.type == "hexagonal" else 4
    side = 2 * dipoles[0]["half_side"]

    out = a.output or os.path.join(a.outdir, f"prismatic_dipoles_{a.count}_{a.type}.yaml")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fmt = lambda v: "[" + " ".join(str(int(x)) for x in v) + "]"
    header = [f"{a.count} {a.type} prismatic dipoles ({nsides} shear loops each) at random positions, "
              f"generated by make_prismatic_dipoles.py",
              f"Loop side {side:.3f} A; separation " +
              (f"{2 * dipoles[0]['half_sep']:.3f} A" if len(a.separation) == 1 else
               f"{min(r['separation'] for r in records):.3f}-{max(r['separation'] for r in records):.3f} A") +
              f"; gap between bounding spheres >= {a.gap} A" + (f"; seed {a.seed}" if a.seed is not None else ""),
              "Per dipole: shear loops, b (a0/2), centre (A), separation (A)"]
    header += [f"  dipole {r['dipole']:>3}: loops {r['shear_loops'][0]}-{r['shear_loops'][-1]}, b {fmt(r['b'])}, "
               f"centre [{' '.join(f'{x:.4f}' for x in r['centre'])}], separation {r['separation']:.3f}"
               for r in records]
    E.write_yaml(out, general, loops, header)
    expected = dict(dipoles=a.count, closed_loops=2 * a.count, burgers_family="1/2<111>", junctions=0,
                    screw_length=0., prismatic_loops=dict(shape="hexagon" if nsides == 6 else "rhombus",
                                                          sides=nsides, side_length=side,
                                                          signs="opposite within each dipole"),
                    per_dipole=records)
    E.write_sidecar(os.path.splitext(out)[0] + ".json", E.load_bcc(), general, loops, expected,
                    extra=dict(placement=dict(type=a.type, cells=cells, gap=a.gap, seed=a.seed,
                                              separation=a.separation, b=a.b, random_sign=a.random_sign),
                               note="shear loops of the same b on the faces of a glide prism; the screw sides "
                                    "cancel pairwise, leaving the two end loops of each dipole"))
    vfrac = sum(4 / 3 * math.pi * r ** 3 for r in radii) / box ** 3
    print(f"{out}: {a.count} {a.type} dipoles, side {side:.3f} A, Loopon_size {da:.6f} A, "
          f"box {cells}^3 cells ({2 * cells ** 3:,} atoms), {len(loops)} shear loops, "
          f"bounding spheres fill {100 * vfrac:.1f} % of the box")
    print("  checks passed: screw sides cancel in every dipole, faces on mid-planes, no overlap")


if __name__ == "__main__":
    main()
