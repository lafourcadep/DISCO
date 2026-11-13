#!/usr/bin/env python3
"""Write a DISCO input with prismatic loop dipoles built from shear loops (BCC, 1/2<111>).

Construction (one dipole). Shear loops with the SAME Burgers vector b lie on the faces of a
glide prism parallel to b:
  rhombic    4 rectangular loops on the two {110} planes n1, n2 containing b (rhombic tube);
  hexagonal  6 rectangular loops on the three {110} planes containing b (hexagonal tube, the
             glide prism of a BCC prismatic loop; the end loops have <112> sides).
Every face normal points out of the tube, so the tube is one consistently oriented surface:
the screw sides (parallel to b) of neighbouring faces coincide with opposite line senses and
cancel. What remains are the edge sides at the two tube ends, normal to b: two prismatic loops
of opposite sign (a prismatic dipole), separated by --separation along b.

Exact cancellation needs the shared sides to coincide, so the requested sizes are adjusted and
the values actually used are printed:
  * the face planes are put on mid-planes between atomic planes, n.x = +/-(k + 1/2) d, so
    DISCO's mid-plane snap does not move them (this sets the side length);
  * Loopon_size is tuned (close to --loopon-size) so that the half-side and half the separation
    are (k + 1/2) Loopon_size: the exact contour with Method: line, and the effective extent
    of the loopon grid with Method: loopons, so both methods cancel the screw sides.

No atom is added or removed: the two loops of a dipole share one glide prism and annihilate by
glide at finite temperature. For prismatic loops that survive MD use make_prismatic_loops.py.

Number and placement:
  --count 1 (default)  one dipole centred at the origin, b = --b (default [1 -1 1]), planes
                       --n1 --n2 (rhombic); box --cells (default 50).
  --count N > 1        N dipoles at random (--seed). b: one of the four <111>, dealt in shuffled
                       rounds of four so the families stay balanced (--b-draw random for an
                       independent draw, --b for one direction), with a random sign
                       (--no-random-sign to keep it); rhombic planes: one of the three {110}
                       pairs containing b; separation: one value or uniform in MIN MAX; position:
                       a random BCC lattice site (cube corner or body centre), so every face stays
                       on its mid-plane. A dipole is redrawn when its bounding sphere comes closer
                       than --gap to another one (periodic minimum image). With --margin, every
                       shear loop stays that far from the box faces (no dipole crosses the
                       periodic boundary; e.g. before opening a free surface for indentation).
                       Box: --cells N | NX NY NZ, --size LX LY LZ (A), or by default the cube in
                       which the bounding spheres fill --fill of the volume.

One sample holds one --type and one --side: the loopon step that makes the screw sides cancel
is set by the side, and no step fits a rhombic and a hexagonal dipole at once (their half-sides
are odd multiples of d/sqrt(3) and of d/(2 sqrt(3)), d = {110} spacing).

Usage:
  python3 make_prismatic_dipoles.py                                   # one hexagonal dipole
  python3 make_prismatic_dipoles.py --type rhombic --side 40 --separation 120 --cells 60
  python3 make_prismatic_dipoles.py --type rhombic --b 1 1 1 --n1 1 -1 0 --n2 0 1 -1
  python3 make_prismatic_dipoles.py --count 10                        # 10 hexagonal dipoles, random b
  python3 make_prismatic_dipoles.py --count 20 --type rhombic --side 40 --separation 60 120 --seed 3
  python3 make_prismatic_dipoles.py --count 12 --side 40 --separation 50 --size 600 600 300 --margin 10 --gap 20
"""
import argparse
import math
import os
import sys

import numpy as np

import exec_input as E

SCRIPT = "make_prismatic_dipoles.py"
DEFAULT_SIDE = {"hexagonal": 60.67, "rhombic": 56.6}


# =============================================================================
# One dipole, centred on the origin
# =============================================================================
def _grid(d, half_side, separation, loopon_target):
    """Loopon step da with half_side = (k_p + 1/2) da, and half the separation (k_m + 1/2) da."""
    k_p = max(1, round(half_side / loopon_target - 0.5))
    da = half_side / (k_p + 0.5)
    k_m = max(1, round(separation / 2 / da - 0.5))
    return da, k_p, k_m, (k_m + 0.5) * da


def design_rhombic(a0, b, n1, n2, side, separation, loopon_target):
    """Four faces +/-n1, +/-n2 of the rhombic tube around b."""
    b, n1, n2 = (np.asarray(v, int) for v in (b, n1, n2))
    for n in (n1, n2):
        if b @ n != 0:
            sys.exit(f"b = {b} does not lie in plane {n} (b.n = {b @ n})")
    if np.linalg.norm(np.cross(n1, n2)) < 1e-12:
        sys.exit("n1 and n2 must be different planes")
    d = E.mid_spacing(n1, a0)
    if abs(d - E.mid_spacing(n2, a0)) > 1e-12:
        sys.exit("n1 and n2 must have the same plane spacing (e.g. both {110})")
    bh = b / np.linalg.norm(b)
    nh1, nh2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
    sin = abs(nh2 @ np.cross(bh, nh1))           # rhombus: half-side = s / sin
    k_s = max(0, round(side / 2 * sin / d - 0.5))
    s = (k_s + 0.5) * d                          # face planes on mid-planes
    half_side = s / sin
    da, k_p, k_m, half_sep = _grid(d, half_side, separation, loopon_target)
    faces = []                                   # (integer outward normal, centre)
    for n_int, other in ((n1, nh2), (n2, nh1)):
        for sgn in (1, -1):
            nh = sgn * n_int / np.linalg.norm(n_int)
            faces.append((sgn * n_int, np.linalg.solve([nh, other, bh], [s, 0., 0.])))  # mid-side, b.x = 0
    return dict(type="rhombic", b=b, bh=bh, d=d, s=s, k_s=k_s, da=da, k_p=k_p, k_m=k_m,
                half_side=half_side, half_sep=half_sep, faces=faces)


def design_hexagonal(a0, b, side, separation, loopon_target):
    """Six faces +/-n1, +/-n2, +/-n3 of the hexagonal tube around b (the three {110} of b)."""
    b = np.asarray(b, int)
    if sorted(np.abs(b)) != [1, 1, 1]:
        sys.exit(f"b = {b}: give a <111> direction (b = a0/2 [b])")
    normals = E.planes_110_containing(b)
    bh = b / np.linalg.norm(b)
    d = E.mid_spacing(normals[0], a0)            # all {110}: same spacing a0/sqrt(2)
    # regular hexagon: apothem s, half-side s / sqrt(3), side 2 s / sqrt(3)
    k_s = max(0, round(side * np.sqrt(3) / 2 / d - 0.5))
    s = (k_s + 0.5) * d
    half_side = s / np.sqrt(3)
    da, k_p, k_m, half_sep = _grid(d, half_side, separation, loopon_target)
    faces = [(sgn * n, s * sgn * n / np.linalg.norm(n)) for n in normals for sgn in (1, -1)]
    return dict(type="hexagonal", b=b, bh=bh, d=d, s=s, k_s=k_s, da=da, k_p=k_p, k_m=k_m,
                half_side=half_side, half_sep=half_sep, faces=faces)


def face_loops(D):
    """The shear loops of one centred dipole, as exec_input loop dicts."""
    return [dict(shape="rectangular", half_m=(D["k_m"] + 0.5) * D["da"], half_p=(D["k_p"] + 0.5) * D["da"],
                 corner_radius=0., center=np.asarray(c, float), m=[int(x) for x in D["b"]],
                 n=[int(x) for x in n_int]) for n_int, c in D["faces"]]


def check_cancellation(D):
    """The screw sides of neighbouring faces coincide with opposite line senses; the
    remaining (edge) sides are normal to b."""
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
    if not (len(screw) == len(edge) == 2 * len(D["faces"])) or not all(
            any(np.allclose(p, b2) and np.allclose(q, a2) for a2, b2 in screw) for p, q in screw):
        sys.exit("internal error: the screw sides do not cancel")


def bounding_radius(D):
    """Largest distance from the dipole centre to a corner of its shear loops."""
    r = 0.
    for n_int, c in D["faces"]:
        p = np.cross(D["bh"], n_int / np.linalg.norm(n_int))
        for u in (1, -1):
            for v in (1, -1):
                r = max(r, np.linalg.norm(c + u * D["half_sep"] * D["bh"] + v * D["half_side"] * p))
    return r


def reach(D):
    """Per-axis half-extent of one centred dipole: how far its shear loops reach from its centre (A)."""
    return np.max([np.abs(l["center"]) + E.loop_half_extent(l) for l in face_loops(D)], axis=0)


def b_directions(a, rng):
    """The Burgers direction of each dipole, before the random sign."""
    if a.b is not None:
        return [np.array(a.b)] * a.count
    if a.b_draw == "random":
        return [E.D111[int(rng.integers(4))] for _ in range(a.count)]
    order = []                                                # shuffled rounds of the four <111>
    while len(order) < a.count:
        r = list(rng.permutation(4))
        if order and r[0] == order[-1]:                       # no repeat across two rounds
            r = r[1:] + r[:1]
        order += r
    return [E.D111[k] for k in order[:a.count]]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--count", type=int, default=1, help="number of dipoles (default 1: one centred dipole)")
    p.add_argument("--type", default="hexagonal", choices=["hexagonal", "rhombic"])
    p.add_argument("--side", type=float, help="side length of the prismatic loops (A; default 60.67 hexagonal, "
                   "56.6 rhombic); adjusted to mid-planes, the same for every dipole")
    p.add_argument("--separation", type=float, nargs="+", default=[80.3], metavar="A",
                   help="distance between the two prismatic loops along b (A): one value, or MIN MAX for a "
                        "uniform random value per dipole; adjusted to the loopon grid")
    p.add_argument("--b", type=int, nargs=3, help="Burgers direction <111> (b = a0/2 [b]) for every dipole "
                   "(default: [1 -1 1] with --count 1, else drawn among the four <111>, see --b-draw)")
    p.add_argument("--n1", type=int, nargs=3, help="--count 1, rhombic: first {110} plane containing b "
                   "(default: 1 1 0 and 0 1 1 for the default b, else two {110} planes of b)")
    p.add_argument("--n2", type=int, nargs=3, help="--count 1, rhombic: second {110} plane containing b")
    p.add_argument("--b-draw", default="balanced", choices=["balanced", "random"],
                   help="balanced: the four <111> dealt in shuffled rounds (equal counts within one); "
                        "random: independent uniform draw per dipole")
    p.add_argument("--random-sign", action=argparse.BooleanOptionalAction, default=True,
                   help="--count > 1: flip the sign of b at random (which end loop is positive); default on")
    box = p.add_mutually_exclusive_group()
    box.add_argument("--cells", type=int, nargs="+", metavar="N",
                     help="box in unit cells: N (cube) or NX NY NZ; default 50 with --count 1, else a cube "
                          "where the dipoles' bounding spheres fill --fill of the volume")
    box.add_argument("--size", type=float, nargs=3, metavar=("LX", "LY", "LZ"),
                     help="box size in A, rounded to whole unit cells")
    p.add_argument("--fill", type=float, default=0.15,
                   help="default box with --count > 1: volume fraction of the bounding spheres (with --gap)")
    p.add_argument("--gap", type=float, default=10.0, help="minimum distance between the bounding spheres "
                   "of two dipoles (A, periodic minimum image)")
    p.add_argument("--margin", type=float,
                   help="keep every shear loop at least this far from the box faces (A), so no dipole "
                        "crosses the periodic boundary; default: no limit")
    p.add_argument("--loopon-size", type=float, default=1.32,
                   help="target loopon grid step (A); tuned slightly so the edges meet exactly")
    p.add_argument("--tries", type=int, default=10000, help="placement attempts per dipole")
    E.add_common_args(p, E.A0_TA, outdir=True)
    a = p.parse_args()

    if a.count < 1:
        sys.exit("--count must be >= 1")
    single = a.count == 1
    if len(a.separation) not in (1, 2) or min(a.separation) <= 0:
        sys.exit("--separation: one positive value, or MIN MAX")
    if len(a.separation) == 2 and a.separation[0] > a.separation[1]:
        sys.exit("--separation MIN MAX: MIN > MAX")
    if a.b is not None and sorted(abs(x) for x in a.b) != [1, 1, 1]:
        sys.exit(f"--b {a.b}: give a <111> direction")
    if (a.n1 is None) != (a.n2 is None):
        sys.exit("--n1 and --n2 go together")
    if a.n1 is not None and (not single or a.type != "rhombic"):
        sys.exit("--n1 --n2: only for one rhombic dipole (--count 1 --type rhombic)")
    cells = E.parse_cells(a.cells)
    if a.margin is not None and a.margin < 0:
        sys.exit("--margin must be >= 0")
    if a.side is None:
        a.side = DEFAULT_SIDE[a.type]
    if single and a.b is None:
        a.b = [1, -1, 1]
    rng = np.random.default_rng(a.seed)

    # ---- the dipoles (b, planes, separation), centred on the origin ----
    dipoles = []
    for b in b_directions(a, rng):
        if not single and a.random_sign and rng.random() < 0.5:
            b = -b
        pair = None
        if a.type == "rhombic":
            normals = E.planes_110_containing(b)
            if a.n1 is not None:
                pair = (a.n1, a.n2)
            elif single:            # (110) / (011) for the default b, as the original rhombic script
                pair = ([1, 1, 0], [0, 1, 1]) if list(b) == [1, -1, 1] else (normals[0], normals[1])
            else:
                i, j = sorted(rng.choice(3, size=2, replace=False))
                pair = (normals[i], normals[j])
        sep = a.separation[0] if len(a.separation) == 1 else float(rng.uniform(*a.separation))
        if a.type == "hexagonal":
            D = design_hexagonal(a.a0, b, a.side, sep, a.loopon_size)
        else:
            D = design_rhombic(a.a0, b, pair[0], pair[1], a.side, sep, a.loopon_size)
        check_cancellation(D)
        dipoles.append(D)
    da = dipoles[0]["da"]
    if len({round(float(D["da"]), 10) for D in dipoles}) > 1:
        sys.exit("internal error: dipoles with different Loopon_size")
    radii = [bounding_radius(D) for D in dipoles]
    reaches = [reach(D) for D in dipoles]
    margin = a.margin or 0.

    # ---- box ----
    if cells is None and a.size:
        cells = tuple(max(1, round(x / a.a0)) for x in a.size)
    elif cells is None and single:
        cells = (50,) * 3
    elif cells is None:
        vol = sum(4 / 3 * math.pi * (r + a.gap / 2) ** 3 for r in radii) / a.fill
        n = math.ceil(vol ** (1 / 3) / a.a0)
        n = max(n, math.ceil(2 * (max(radii) + a.gap) / a.a0))
        if a.margin is not None:
            n = max(n, math.ceil(2 * (np.max(reaches) + margin) / a.a0) + 1)
        cells = (n,) * 3
    box = a.a0 * np.array(cells, float)
    o = -box / 2                                              # DISCO's lattice origin
    if single:
        if np.any(reaches[0] + max(margin, 10.) > box / 2):
            sys.exit(f"the dipole reaches {np.round(reaches[0], 1)} A from the centre, the half-box is "
                     f"{np.round(box / 2, 1)} A: increase --cells (keeping {max(margin, 10.)} A free)")
    else:
        if 2 * (max(radii) + a.gap / 2) > box.min():
            sys.exit(f"a dipole spans {2 * max(radii):.1f} A (+ --gap {a.gap} A), more than the {box.min():.1f} A "
                     f"box: increase --cells or --size")
        if a.margin is not None and np.any(np.max(reaches, axis=0) + margin >= box / 2):
            sys.exit(f"a dipole reaches {np.round(np.max(reaches, axis=0), 1)} A from its centre; with "
                     f"--margin {margin} it does not fit in the {np.round(box, 1)} A box")

    # ---- placement: on lattice vectors, no overlap of bounding spheres ----
    centres = [np.zeros(3)] if single else [None] * a.count
    for t in ([] if single else np.argsort(radii)[::-1]):       # large ones first
        lo, hi = (o + reaches[t] + margin, -o - reaches[t] - margin) if a.margin is not None else (o, -o)
        for _ in range(a.tries):
            c = o + a.a0 * np.round((rng.uniform(lo, hi) - o) / a.a0)   # cube-corner site
            if rng.random() < 0.5:
                c = c + 0.5 * a.a0                                     # or body-centre site
            if a.margin is not None and (np.any(c < lo - 1e-9) or np.any(c > hi + 1e-9)):
                continue
            if all(np.linalg.norm(E.min_image(c - centres[u], box)) >= radii[t] + radii[u] + a.gap
                   for u in range(a.count) if centres[u] is not None):
                centres[t] = c
                break
        else:
            placed = sum(c is not None for c in centres)
            sys.exit(f"placed only {placed} of {a.count} dipoles ({a.tries} tries for the next one): increase "
                     f"--cells, lower --fill, --gap or --margin, or reduce --count")

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
    seps = [r["separation"] for r in records]
    header = [f"{a.count} {a.type} prismatic dipole(s) ({nsides} shear loops each), generated by {SCRIPT}",
              f"Loop side {side:.3f} A; separation "
              + (f"{seps[0]:.3f} A" if len(a.separation) == 1 else f"{min(seps):.3f}-{max(seps):.3f} A")
              + ("; centred" if single else f"; random positions, gap between bounding spheres >= {a.gap} A"
                 + (f"; >= {margin} A from the box faces" if a.margin is not None else "") + f"; seed {a.seed}"),
              "Per dipole: shear loops, b (a0/2), centre (A), separation (A)"]
    header += [f"  dipole {r['dipole']:>3}: loops {r['shear_loops'][0]}-{r['shear_loops'][-1]}, b {fmt(r['b'])}, "
               f"centre [{' '.join(f'{x:.4f}' for x in r['centre'])}], separation {r['separation']:.3f}"
               for r in records]
    E.write_yaml(out, general, loops, header)
    expected = dict(closed_loops=2 * a.count, burgers_families=["1/2<111>"], junctions=0, dipoles=a.count,
                    screw_length=0., prismatic_loops=dict(shape="hexagon" if nsides == 6 else "rhombus",
                                                          sides=nsides, side_length=side,
                                                          signs="opposite within each dipole"),
                    per_dipole=records)
    E.write_sidecar(os.path.splitext(out)[0] + ".json", E.load_bcc(), general, loops, expected, script=SCRIPT,
                    extra=dict(placement=dict(type=a.type, cells=list(cells), gap=a.gap, margin=a.margin,
                                              seed=None if single else a.seed, separation=a.separation, b=a.b,
                                              b_draw=None if a.b else a.b_draw,
                                              random_sign=a.random_sign and not single),
                               note="shear loops of the same b on the faces of a glide prism; the screw sides "
                                    "cancel pairwise, leaving the two end loops of each dipole"))
    vfrac = sum(4 / 3 * math.pi * r ** 3 for r in radii) / box.prod()
    D0 = dipoles[0]
    print(f"{out}: {a.count} {a.type} dipole(s), {len(loops)} shear loops, "
          f"box {' x '.join(map(str, cells))} cells ({2 * int(np.prod(cells)):,} atoms)")
    print(f"  loop side      {side:8.3f} A  (asked {a.side}; face planes at n.x = +/-{D0['s']:.3f} A "
          f"= ({D0['k_s']}+1/2) d, d = {D0['d']:.4f} A)")
    print(f"  separation     " + (f"{seps[0]:8.3f} A" if len(seps) == 1 else f"{min(seps):.3f}-{max(seps):.3f} A")
          + f"  (asked {' '.join(map(str, a.separation))})")
    print(f"  Loopon_size    {da:8.6f} A  (asked {a.loopon_size})")
    if not single:
        fam = {fmt(v): sum(int(np.array_equal(np.array(r["b"]) * np.sign(r["b"][0]), v)) for r in records)
               for v in E.D111}
        print(f"  b families {fam}; bounding spheres fill {100 * vfrac:.1f} % of the box")
    print("  checks passed: screw sides cancel in every dipole, faces on mid-planes"
          + ("" if single else ", no overlap") + (", loops inside the margin" if a.margin is not None else ""))


if __name__ == "__main__":
    main()
