#!/usr/bin/env python3
"""Write an exec_cuda input with N prismatic loop dipoles placed at random in the box.

Each dipole is built as in make_prismatic_dipole.py (rhombic: 4 shear loops on the two
{110} planes of one pair containing b) or make_hexagonal_prismatic_dipole.py (hexagonal:
6 shear loops on the three {110} planes containing b): the screw sides cancel pairwise and
leave two prismatic loops of opposite sign, normal to b, separated along b.

Per dipole, at random (--seed):
  b             one of the four <111> directions, dealt in shuffled rounds of four so the
                families stay balanced (12 dipoles = 3 per family; --b-draw random for an
                independent draw, or --b for one direction), with a random sign
                (--random-sign), so the loop of positive sign is at +b or -b;
  planes        rhombic: one of the three {110} pairs containing b;
  separation    uniform in [--separation MIN MAX], on the loopon grid ((k + 1/2) Loopon_size);
  position      a random BCC lattice site of the box (cube corner or body centre), so
                every face plane stays on its mid-plane and exec_cuda's snap does not move it.
A dipole is rejected and redrawn when its bounding sphere comes closer than --gap to that of
another dipole (periodic minimum image). With --margin, every shear loop also stays at least
that far from the box faces, so no dipole crosses the periodic boundary (needed when a free
surface is opened later, e.g. for indentation).

The box is a cube by default (--cells N, or sized from --fill); --cells NX NY NZ or
--size LX LY LZ give any orthorhombic box.

One sample holds one type: the loopon step that makes the screw sides cancel exactly is
set by the side length, and no step fits a rhombic and a hexagonal dipole at once (their
half-sides are odd multiples of d/sqrt(3) and of d/(2 sqrt(3)), d = {110} spacing). All
dipoles of a sample therefore have the same side; the separation can differ.

Usage:
  python3 make_prismatic_dipoles.py --count 10                       # 10 hexagonal dipoles, random b
  python3 make_prismatic_dipoles.py --count 20 --type rhombic --side 40 --separation 60 120 --seed 3
  python3 make_prismatic_dipoles.py --count 8 --b 1 1 1 --cells 120 --format lammps
  python3 make_prismatic_dipoles.py --count 12 --side 40 --separation 50 --size 600 600 300 --margin 10 --gap 20
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


def reach(D):
    """Per-axis half-extent of one centred dipole: how far its shear loops reach from its centre (A)."""
    return np.max([np.abs(l["center"]) + E.loop_half_extent(l) for l in face_loops(D)], axis=0)


def b_directions(a, rng):
    """The Burgers direction of each dipole, before the random sign."""
    if a.b is not None:
        return [np.array(a.b)] * a.count
    four = [v for v in B111 if v[0] > 0]
    if a.b_draw == "random":
        return [four[int(rng.integers(4))] for _ in range(a.count)]
    order = []                                                # shuffled rounds of the four <111>
    while len(order) < a.count:
        r = list(rng.permutation(4))
        if order and r[0] == order[-1]:                       # no repeat across two rounds
            r = r[1:] + r[:1]
        order += r
    return [four[k] for k in order[:a.count]]


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
                   "(default: drawn among the four <111>, see --b-draw)")
    p.add_argument("--b-draw", default="balanced", choices=["balanced", "random"],
                   help="balanced: the four <111> dealt in shuffled rounds (equal counts within one); "
                        "random: independent uniform draw per dipole")
    p.add_argument("--random-sign", action=argparse.BooleanOptionalAction, default=True,
                   help="flip the sign of b at random (which end loop is positive); default on")
    box = p.add_mutually_exclusive_group()
    box.add_argument("--cells", type=int, nargs="+", metavar="N",
                     help="box in unit cells: N (cube) or NX NY NZ; default: a cube where the dipoles' "
                          "bounding spheres fill --fill of the volume")
    box.add_argument("--size", type=float, nargs=3, metavar=("LX", "LY", "LZ"),
                     help="box size in A, rounded to whole unit cells")
    p.add_argument("--fill", type=float, default=0.15,
                   help="default box: volume fraction of the bounding spheres (with --gap)")
    p.add_argument("--gap", type=float, default=10.0, help="minimum distance between the bounding spheres "
                   "of two dipoles (A, periodic minimum image)")
    p.add_argument("--margin", type=float,
                   help="keep every shear loop at least this far from the box faces (A), so no dipole "
                        "crosses the periodic boundary; default: no limit")
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
    if a.cells is not None and (len(a.cells) not in (1, 3) or min(a.cells) < 1):
        sys.exit("--cells: N, or NX NY NZ, positive")
    if a.margin is not None and a.margin < 0:
        sys.exit("--margin must be >= 0")
    if a.side is None:
        a.side = 60.67 if a.type == "hexagonal" else 56.6
    rng = np.random.default_rng(a.seed)

    # ---- draw the dipoles (b, planes, separation), centred on the origin ----
    dipoles = []
    for b in b_directions(a, rng):
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
    reaches = [reach(D) for D in dipoles]
    margin = a.margin or 0.

    # ---- box ----
    if a.cells:
        cells = tuple(a.cells) * 3 if len(a.cells) == 1 else tuple(a.cells)
    elif a.size:
        cells = tuple(max(1, round(x / a.a0)) for x in a.size)
    else:
        vol = sum(4 / 3 * math.pi * (r + a.gap / 2) ** 3 for r in radii) / a.fill
        n = math.ceil(vol ** (1 / 3) / a.a0)
        n = max(n, math.ceil(2 * (max(radii) + a.gap) / a.a0))
        if a.margin is not None:
            n = max(n, math.ceil(2 * (np.max(reaches) + margin) / a.a0) + 1)
        cells = (n,) * 3
    box = a.a0 * np.array(cells, float)
    o = -box / 2                                              # exec_cuda's lattice origin
    if 2 * (max(radii) + a.gap / 2) > box.min():
        sys.exit(f"a dipole spans {2 * max(radii):.1f} A (+ --gap {a.gap} A), more than the {box.min():.1f} A box: "
                 f"increase --cells or --size")
    if a.margin is not None and np.any(np.max(reaches, axis=0) + margin >= box / 2):
        sys.exit(f"a dipole reaches {np.round(np.max(reaches, axis=0), 1)} A from its centre; with "
                 f"--margin {margin} it does not fit in the {np.round(box, 1)} A box")

    # ---- random placement on lattice vectors, no overlap of bounding spheres ----
    order = np.argsort(radii)[::-1]                           # large ones first
    centres = [None] * a.count
    for t in order:
        lo, hi = (o + reaches[t] + margin, -o - reaches[t] - margin) if a.margin is not None else (o, -o)
        for _ in range(a.tries):
            c = o + a.a0 * np.round((rng.uniform(lo, hi) - o) / a.a0)   # cube-corner site
            if rng.random() < 0.5:
                c = c + 0.5 * a.a0                                     # or body-centre site
            if a.margin is not None and (np.any(c < lo - 1e-9) or np.any(c > hi + 1e-9)):
                continue
            if all(np.linalg.norm(min_image(c - centres[u], box)) >= radii[t] + radii[u] + a.gap
                   for u in range(a.count) if centres[u] is not None):
                centres[t] = c
                break
        else:
            placed = sum(c is not None for c in centres)
            sys.exit(f"placed only {placed} of {a.count} dipoles ({a.tries} tries for the next one): increase --cells, "
                     f"lower --fill, --gap or --margin, or reduce --count")

    # ---- loops, checks, records ----
    loops, records = [], []
    for t, (D, c) in enumerate(zip(dipoles, centres)):
        first = len(loops) + 1
        for l in face_loops(D):
            l["center"] = l["center"] + c
            moved = np.linalg.norm(E.snap_center(l["center"], l["n"], a.a0, cells) - l["center"])
            if moved > 1e-6:
                sys.exit(f"internal error: dipole {t + 1}: a face is not on a mid-plane (moved {moved:.2e} A)")
            if a.margin is not None and np.any(np.abs(l["center"]) + E.loop_half_extent(l) > box / 2 - margin + 1e-9):
                sys.exit(f"internal error: dipole {t + 1}: a loop is closer than --margin to the box faces")
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
              f"; gap between bounding spheres >= {a.gap} A" +
              (f"; >= {margin} A from the box faces" if a.margin is not None else "") + (f"; seed {a.seed}" if a.seed is not None else ""),
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
                    extra=dict(placement=dict(type=a.type, cells=list(cells), gap=a.gap, margin=a.margin,
                                              seed=a.seed, separation=a.separation, b=a.b,
                                              b_draw=None if a.b else a.b_draw, random_sign=a.random_sign),
                               note="shear loops of the same b on the faces of a glide prism; the screw sides "
                                    "cancel pairwise, leaving the two end loops of each dipole"))
    vfrac = sum(4 / 3 * math.pi * r ** 3 for r in radii) / box.prod()
    fam = {fmt(v): sum(int(np.array_equal(np.array(r["b"]) * np.sign(r["b"][0]), v)) for r in records) for v in B111 if v[0] > 0}
    print(f"{out}: {a.count} {a.type} dipoles, side {side:.3f} A, Loopon_size {da:.6f} A, "
          f"box {' x '.join(map(str, cells))} cells ({2 * int(np.prod(cells)):,} atoms), {len(loops)} shear loops, "
          f"bounding spheres fill {100 * vfrac:.1f} % of the box")
    print(f"  b families {fam}")
    print("  checks passed: screw sides cancel in every dipole, faces on mid-planes, no overlap" +
          (", loops inside the margin" if a.margin is not None else ""))


if __name__ == "__main__":
    main()
