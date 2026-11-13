#!/usr/bin/env python3
"""Write a DISCO input with N FCC dislocation loops made by adding or removing atoms:
perfect prismatic loops on their parallelepiped glide prism, Frank loops, or a mix of both
(an irradiation-like population of interstitial and vacancy loops).

Loop kinds (--kind prismatic | frank | mixed, mixed: Frank with probability --frank-fraction):
  prismatic  perfect b = 1/2<110>, loop plane normal to b ({110}): pure-edge prismatic loop,
             2 {110} atomic planes added or removed. Exactly two {111} planes contain b, so
             the loop is a RHOMBUS (default --shape) whose sides are their intersections with
             the {110} plane: <112> directions at 70.5 / 109.5 deg. Every side can glide along
             b on its own {111} face, so the loop moves rigidly along a parallelepiped prism
             (the FCC counterpart of the hexagonal 1/2<111> prism in BCC). --shape circular
             gives a round loop on the same plane instead.
  frank      faulted b = a0/3<111>, loop on the {111} plane normal to b: sessile, one {111}
             layer added (extrinsic fault) or removed (intrinsic fault). --shape hexagonal
             (default) or triangular: sides along <110>, the intersections with the other
             {111} planes, as in quenched and irradiated FCC metals; or circular.
Per loop:
  b       random direction (6 <110> or 4 <111>, random sign), or --b (one direction for all
          loops of that kind; a <110> for prismatic, a <111> for Frank);
  type    interstitial (m.n < 0) or vacancy (m.n > 0): --type interstitial|vacancy|mixed
          (mixed: random per loop);
  size    --size: rhombus side, hexagon / triangle circumradius, or circle radius (A);
  centre  random in the box (--seed). A loop is redrawn if its bounding sphere comes closer
          than --gap to another one (periodic minimum image) and, for two perfect loops of
          OPPOSITE type with the same b direction, if their glide prisms overlap (lateral
          distance perpendicular to b < sum of the bounding radii + --gap): such a pair would
          glide together and annihilate. Frank loops are sessile and only need --gap.

The YAML uses Lattice: FCC, Element (--element, default Cu; Mass from DISCO's table or --mass)
and Method: line. Burgers vectors of Frank loops are written as exact fractions (2/3).
The JSON sidecar holds the ground truth: Burgers vector, kind, type, glide planes and the
expected number of atoms added / removed per loop. There is no DANCE slip-system table for FCC,
so expected_burgers_type / expected_slip_system are null.

Usage:
  python3 make_fcc_loops.py                                          # 1 interstitial rhombus, side 25 A
  python3 make_fcc_loops.py --count 8 --type mixed --seed 2          # 8 prismatic rhombi
  python3 make_fcc_loops.py --kind frank --shape triangular --type vacancy --size 20
  python3 make_fcc_loops.py --kind mixed --count 12 --type mixed --frank-fraction 0.5 --cells 60
  python3 make_fcc_loops.py --kind prismatic --b 1 -1 0 --size 30 --element Ni --a0 3.52
"""
import argparse
import math
import os
import sys

import numpy as np

import exec_input as E

SCRIPT = "make_fcc_loops.py"
SHAPES = {"prismatic": ("rhombus", "circular"), "frank": ("hexagonal", "triangular", "circular")}


def triangle_vertex_dir(n):
    """A <112> in the {111} plane n: the triangle vertices along <112> put its sides on <110>."""
    u = next(v for v in E.D110 if v @ n == 0)
    return [int(x) for x in E.reduced(np.cross(n, u))]


def shape_of(kind, shape, size, b):
    if shape == "circular":
        return dict(shape="circular", radius=size)
    if kind == "prismatic":
        return E.rhombus_110(b, size)
    if shape == "hexagonal":                      # default in-plane axis: a <110> vertex,
        return dict(shape="hexagonal", circumradius=size)   # so the sides run along <110>
    return dict(shape="polygon", sides=3, circumradius=size, vertex_dir=triangle_vertex_dir(b))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--kind", default="prismatic", choices=["prismatic", "frank", "mixed"])
    p.add_argument("--frank-fraction", type=float, default=0.5, help="--kind mixed: probability of a Frank loop")
    p.add_argument("--count", type=int, default=1)
    p.add_argument("--type", default="interstitial", choices=["interstitial", "vacancy", "mixed"])
    p.add_argument("--shape", default=None,
                   help="prismatic: rhombus (default) | circular; frank: hexagonal (default) | triangular | "
                        "circular. With --kind mixed it applies to the kind that accepts it")
    p.add_argument("--size", type=float, default=25.0,
                   help="rhombus side, hexagon / triangle circumradius or circle radius (A)")
    p.add_argument("--b", type=int, nargs=3, help="Burgers direction: a <110> (prismatic) or a <111> (Frank)")
    p.add_argument("--cells", type=int, nargs="+", default=None,
                   help="Nx [Ny Nz] (default: smallest even cube holding the loops)")
    p.add_argument("--gap", type=float, default=15.0, help="minimum distance between loops (A)")
    p.add_argument("--element", default="Cu", help="species name in the output files (default Cu)")
    p.add_argument("--mass", type=float, default=None, help="Dalton (default: DISCO's table for --element)")
    E.add_common_args(p, E.A0_CU, output="fcc_loops.yaml")
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)

    kinds = ["prismatic", "frank"] if a.kind == "mixed" else [a.kind]
    shape_for = {}
    for k in kinds:
        if a.shape is None or (a.kind == "mixed" and a.shape not in SHAPES[k]):
            shape_for[k] = SHAPES[k][0]
        elif a.shape in SHAPES[k]:
            shape_for[k] = a.shape
        else:
            sys.exit(f"--shape {a.shape} is not available for {k} loops: {', '.join(SHAPES[k])}")
    if a.b is not None:
        bb = sorted(abs(x) for x in a.b)
        if bb not in ([0, 1, 1], [1, 1, 1]):
            sys.exit("--b must be a <110> or a <111> direction")
        if a.kind != "mixed" and (bb == [1, 1, 1]) != (a.kind == "frank"):
            sys.exit(f"--b {a.b}: {a.kind} loops need a " + ("<111>" if a.kind == "frank" else "<110>"))
    if not 0. <= a.frank_fraction <= 1.:
        sys.exit("--frank-fraction must be in [0, 1]")

    reach_max = max(E.bounding_radius(shape_of(k, shape_for[k], a.size,
                                               E.D110[0] if k == "prismatic" else E.D111[0])) for k in kinds)
    cells = E.parse_cells(a.cells)
    if cells is None:
        side = (2 * reach_max + a.gap) * max(1., a.count ** (1 / 3)) * 1.6 + 6 * a.a0
        n = math.ceil(side / a.a0)
        cells = [n + n % 2] * 3
    box = a.a0 * np.asarray(cells, float)

    loops, placed = [], []                       # placed: (centre, b_hat, type, kind, reach)
    for k in range(a.count):
        for attempt in range(10000):
            kind = kinds[0] if len(kinds) == 1 else ("frank" if rng.random() < a.frank_fraction else "prismatic")
            family = E.D110 if kind == "prismatic" else E.D111
            if a.b is not None and sorted(abs(x) for x in a.b) == sorted(abs(int(x)) for x in family[0]):
                b = np.array(a.b)
            else:
                b = family[rng.integers(len(family))] * rng.choice((1, -1))
            typ = a.type if a.type != "mixed" else ("interstitial", "vacancy")[rng.integers(2)]
            shape = shape_of(kind, shape_for[kind], a.size, b)
            reach = E.bounding_radius(shape)
            c = (rng.random(3) - 0.5) * box if a.count > 1 or a.cells else np.zeros(3)
            bh = b / np.linalg.norm(b)
            ok = True
            for c2, bh2, typ2, kind2, reach2 in placed:
                r = E.min_image(c - c2, box)
                if np.linalg.norm(r) < reach + reach2 + a.gap:
                    ok = False
                elif kind == kind2 == "prismatic" and typ2 != typ and abs(abs(bh @ bh2) - 1) < 1e-9 and \
                        np.linalg.norm(r - (r @ bh) * bh) < reach + reach2 + a.gap:
                    ok = False                   # same glide prism, opposite sign: would annihilate
                if not ok:
                    break
            if ok:
                break
        else:
            sys.exit(f"could not place loop {k + 1}: increase --cells or reduce --count/--size/--gap")
        placed.append((c, bh, typ, kind, reach))
        # n = b direction; m.n < 0: interstitial. Perfect: m = <110> (b = a0/2 m); Frank: m = 2/3 <111>
        scale = 1. if kind == "prismatic" else 2. / 3.
        m = (-scale if typ == "interstitial" else scale) * b
        loops.append(dict(**shape, center=list(c), m=[float(x) for x in m], n=[int(x) for x in b],
                          kind=kind, type=typ))

    # The box is periodic: a loop may cross a face, but not come within --gap of its own image
    if 2 * reach_max + a.gap > box.min():
        sys.exit(f"box {box.min():.1f} A is too small for loops of reach {reach_max:.1f} A and --gap {a.gap}: "
                 f"increase --cells")
    d_min = min(E.mid_spacing(l["n"], a.a0, "FCC") for l in loops)
    general = E.general_block(a.a0, tuple(cells), 0.33 * d_min, a.format,
                              lattice="FCC", element=a.element, mass=a.mass)
    n_int = sum(l["type"] == "interstitial" for l in loops)
    n_frank = sum(l["kind"] == "frank" for l in loops)
    header = [f"{len(loops)} FCC loop(s) made by adding/removing atoms; generated by {SCRIPT}",
              f"{len(loops) - n_frank} perfect prismatic 1/2<110> ({shape_for.get('prismatic', '-')}), "
              f"{n_frank} Frank a0/3<111> ({shape_for.get('frank', '-')}); size {a.size} A",
              f"{n_int} interstitial, {len(loops) - n_int} vacancy (m.n < 0: interstitial, m.n > 0: vacancy)"]
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    E.write_yaml(a.output, general, loops, header)

    rho = 4 / a.a0 ** 3
    per_loop = []
    for l in loops:
        area = E.contour_area(l)
        nh = np.asarray(l["n"], float) / np.linalg.norm(l["n"])
        height = abs(a.a0 / 2 * np.asarray(l["m"]) @ nh)          # |b.n_hat|
        rec = dict(kind=l["kind"], type=l["type"], normal=l["n"],
                   burgers=("1/2[" if l["kind"] == "prismatic" else "1/3[")      # b = a0/2 m
                   + " ".join(str(int(round(x * (1 if l["kind"] == "prismatic" else 1.5)))) for x in l["m"]) + "]",
                   planes=int(round(height / E.mid_spacing(l["n"], a.a0, "FCC"))),
                   atoms_added=(1 if l["type"] == "interstitial" else -1) * round(rho * area * height))
        if l["kind"] == "prismatic":
            rec["glide_planes"] = [[int(x) for x in q] for q in E.glide_planes_111(np.asarray(l["n"]))]
        else:
            rec["stacking_fault"] = "extrinsic" if l["type"] == "interstitial" else "intrinsic"
        per_loop.append(rec)
    families = sorted({"1/2<110>" if l["kind"] == "prismatic" else "1/3<111>" for l in loops})
    expected = dict(closed_loops=len(loops), burgers_families=families, junctions=0, loops=per_loop)
    E.write_sidecar(os.path.splitext(a.output)[0] + ".json", None, general, loops, expected,
                    extra=dict(note="FCC loops with b.n != 0; DISCO adds/removes the slab atoms. "
                                    "atoms_added is rho A |b.n| (the code's count differs by the "
                                    "atoms cut by the contour)"),
                    script=SCRIPT)
    print(f"wrote {a.output}: {len(loops)} loop(s) ({len(loops) - n_frank} prismatic, {n_frank} Frank; "
          f"{n_int} interstitial), cells {cells}")


if __name__ == "__main__":
    main()
