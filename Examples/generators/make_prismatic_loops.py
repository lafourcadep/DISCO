#!/usr/bin/env python3
"""Write an exec_cuda input with N prismatic dislocation loops made by adding or removing atoms.

Unlike make_prismatic_dipole(s).py (shear loops on the faces of a glide prism, no matter
added or removed, so the two loops of opposite sign glide together and annihilate at finite
temperature), every loop here has its Burgers vector out of its plane (b.n != 0) and
exec_cuda adds (interstitial) or removes (vacancy) the slab of |b.n| atoms inside it
(Code/prismatic.hpp). A single loop, or loops of the same sign, are then stable at
temperature: they can glide along b but cannot annihilate without climb.

Per loop:
  b       a 1/2<111> direction (--b, or random among the four <111>, random sign);
  habit   111: n // b, pure-edge prismatic loop (3 {111} planes in the slab);
          110: n one of the three {110} planes at 35.3 deg from b (mixed loop, 1 plane);
  type    interstitial (m.n < 0) or vacancy (m.n > 0), --type interstitial|vacancy|mixed
          (mixed: random per loop);
  shape   --shape circular (radius --radius, A) or hexagonal (circumradius --radius,
          vertices along --vertex-dir, default the code's in-plane axis: a <110> for {111},
          so the sides run along <112>); --corner-radius rounds the hexagon corners;
  centre  random in the box (--seed). A loop is redrawn if its bounding sphere comes closer
          than --gap to another one (periodic minimum image) and, for two loops of OPPOSITE
          type with the same b direction, if their glide prisms overlap (lateral distance
          perpendicular to b < 2 R + --gap): such a pair would meet by glide and annihilate.

The file uses Method: line (exact field, no loopon limit). Loopon_size, written for a switch to
Method: loopons, defaults to 0.33 d (d: atomic-plane spacing of the habit plane), the largest
step for which the planes next to the cut of an interstitial loop are then smooth.

Usage:
  python3 make_prismatic_loops.py                                   # 1 interstitial loop, R = 20 A
  python3 make_prismatic_loops.py --count 6 --type mixed --radius 15 --cells 50 --seed 2
  python3 make_prismatic_loops.py --b 1 1 1 --type vacancy --habit 110 -o vac110.yaml
  python3 make_prismatic_loops.py --shape hexagonal --radius 22 --corner-radius 2 -o hex.yaml
"""
import argparse
import itertools
import math
import os
import sys

import numpy as np

import exec_input as E

B111 = [np.array(v) for v in itertools.product((1, -1), repeat=3) if v[0] == 1]


def habit_normal(b, habit, rng):
    """Integer loop normal n with n.b > 0: b itself ({111}) or a {110} at 35.3 deg from b."""
    if habit == "111":
        return np.array(b)
    cands = [np.array(n) for n in set(itertools.permutations((1, 1, 0))) | set(itertools.permutations((1, -1, 0)))
             if np.dot(n, b) == 2]
    return cands[rng.integers(len(cands))]


def min_image(r, box):
    return r - box * np.round(r / box)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--output", default="prismatic_loops.yaml")
    p.add_argument("--count", type=int, default=1)
    p.add_argument("--type", default="interstitial", choices=["interstitial", "vacancy", "mixed"])
    p.add_argument("--habit", default="111", choices=["111", "110"])
    p.add_argument("--shape", default="circular", choices=["circular", "hexagonal"])
    p.add_argument("--radius", type=float, default=20.0, help="loop radius, or hexagon circumradius (A)")
    p.add_argument("--corner-radius", type=float, default=0.0, help="hexagon corner rounding (A)")
    p.add_argument("--vertex-dir", type=int, nargs=3, default=None,
                   help="hexagon: centre -> first vertex, in the plane (default: code's in-plane axis)")
    p.add_argument("--b", type=int, nargs=3, help="Burgers direction <111> (default: random per loop)")
    p.add_argument("--cells", type=int, nargs="+", default=None,
                   help="Nx [Ny Nz] (default: smallest even cube holding the loops)")
    p.add_argument("--gap", type=float, default=15.0, help="minimum distance between loops (A)")
    p.add_argument("--loopon-size", type=float, default=None, help="default 0.33 d")
    p.add_argument("--a0", type=float, default=E.A0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--format", default="lammps", choices=["xyz", "lammps", "both", "none"])
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)

    if a.cells is None:
        side = (2 * a.radius + a.gap) * max(1., a.count ** (1 / 3)) * 1.6 + 6 * a.a0
        n = math.ceil(side / a.a0)
        cells = [n + n % 2] * 3
    else:
        cells = (a.cells * 3)[:3] if len(a.cells) == 1 else a.cells
        if len(cells) != 3:
            sys.exit("--cells takes 1 or 3 integers")
    box = a.a0 * np.asarray(cells, float)

    loops, placed = [], []          # placed: (centre, b_hat, sign)
    for k in range(a.count):
        for attempt in range(10000):
            b = np.array(a.b) if a.b else B111[rng.integers(4)] * rng.choice((1, -1))
            if sorted(abs(int(x)) for x in b) != [1, 1, 1]:
                sys.exit("--b must be a <111> direction")
            n = habit_normal(b, a.habit, rng)
            kind = a.type if a.type != "mixed" else ("interstitial", "vacancy")[rng.integers(2)]
            m = -b if kind == "interstitial" else b          # m.n < 0: interstitial
            c = (rng.random(3) - 0.5) * box if a.count > 1 or a.cells else np.zeros(3)
            bh = b / np.linalg.norm(b)
            ok = True
            for c2, bh2, kind2 in placed:
                r = min_image(c - c2, box)
                if np.linalg.norm(r) < 2 * a.radius + a.gap:
                    ok = False
                elif kind2 != kind and abs(abs(bh @ bh2) - 1) < 1e-9 and \
                        np.linalg.norm(r - (r @ bh) * bh) < 2 * a.radius + a.gap:
                    ok = False                                  # same glide prism, opposite sign
                if not ok:
                    break
            if ok:
                break
        else:
            sys.exit(f"could not place loop {k + 1}: increase --cells or reduce --count/--radius/--gap")
        placed.append((c, bh, kind))
        if a.shape == "hexagonal":
            if a.vertex_dir is not None and np.dot(a.vertex_dir, n) != 0:
                sys.exit(f"--vertex-dir {a.vertex_dir} is not in the plane n = {list(n)} of loop {k + 1}")
            shape = dict(shape="hexagonal", circumradius=a.radius, vertex_dir=a.vertex_dir,
                         corner_radius=a.corner_radius)
        else:
            shape = dict(shape="circular", radius=a.radius)
        loops.append(dict(**shape, center=list(c), m=[int(x) for x in m],
                          n=[int(x) for x in n], kind=kind))

    d = E.mid_spacing(loops[0]["n"], a.a0)
    da = a.loopon_size or 0.33 * d
    general = E.general_block(a.a0, tuple(cells), da, a.format)
    n_int = sum(l["kind"] == "interstitial" for l in loops)
    header = [f"{len(loops)} prismatic loop(s) made by adding/removing atoms; generated by make_prismatic_loops.py",
              f"{n_int} interstitial, {len(loops) - n_int} vacancy; habit {{{a.habit}}}, {a.shape}, R = {a.radius} A"
              + (f", corner radius {a.corner_radius} A" if a.shape == "hexagonal" and a.corner_radius else ""),
              "m.n < 0: interstitial (atoms added), m.n > 0: vacancy (atoms removed)"]
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    E.write_yaml(a.output, general, loops, header)

    if a.shape == "hexagonal":    # rounded regular hexagon, r clamped to the apothem
        rc = min(a.corner_radius, a.radius * math.sqrt(3) / 2)
        area = 1.5 * math.sqrt(3) * a.radius ** 2 - rc ** 2 * (2 * math.sqrt(3) - math.pi)
    else:
        area = math.pi * a.radius ** 2
    rho = 2 / a.a0 ** 3
    expected = dict(closed_loops=len(loops), burgers_family="1/2<111>", junctions=0,
                    prismatic_loops=[dict(type=l["kind"], normal=l["n"],
                                          atoms_added=(1 if l["kind"] == "interstitial" else -1)
                                          * round(rho * area * abs(
                                              np.dot(l["m"], l["n"])) * a.a0 / 2 / np.linalg.norm(l["n"])))
                                     for l in loops])
    E.write_sidecar(os.path.splitext(a.output)[0] + ".json", E.load_bcc(), general, loops, expected,
                    extra=dict(note="loops with b.n != 0; exec_cuda adds/removes the slab atoms"),
                    script="make_prismatic_loops.py")
    print(f"wrote {a.output}: {len(loops)} loop(s) ({n_int} interstitial), cells {cells}, "
          f"Loopon_size {da:.4f} A (d = {d:.4f} A)")


if __name__ == "__main__":
    main()
