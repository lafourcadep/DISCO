#!/usr/bin/env python3
"""Write a DISCO input for a periodic Voronoi polycrystal populated with dislocation loops.

The grains are drawn here (seeds uniform in the box, uniformly random orientations) and written as
an explicit Grain_list (Center + Bunge Euler angles), so the loops can be placed knowing every
grain. Loops are added until the nominal line length reaches --density (m^-2) or --loops; each loop:
  centre  uniform in the box; its grain is the nearest seed (periodic), so grains get loops in
          proportion to their volume. Loops are kept --gap apart (bounding spheres, periodic).
  kind    drawn from --kinds (name:weight list):
            FCC  glide      perfect 1/2<110> on {111}, circular
                 shockley   a0/6<112> on {111}, circular (only the stable sense of b for n, which
                            makes an intrinsic stacking fault)
                 frank      a0/3<111> on {111}, hexagonal with <110> sides (interstitial or vacancy)
                 prismatic  perfect 1/2<110> normal to b, rhombus with <112> sides on its two {111}
                            glide planes (interstitial or vacancy)
            BCC  glide      1/2<111> on {110}, circular
                 prismatic  1/2<111> normal to b, hexagon with <112> sides (on its {110} glide prism)
  size    lognormal: median --size (radius, or the circumradius / reach of the polygon), log-sd
          --size-sigma;
  type    --type interstitial | vacancy | mixed for loops with b.n != 0.
DISCO clips every loop to its grain (a loop larger than its grain spans the grain's cross-section)
and skips a loop that misses its grain; it writes the final list to Poly_..._loops.txt. Loops whose
bounding sphere leaves the grain's cell are marked "may_clip" in the sidecar.

Usage:
  python3 make_polycrystal_loops.py                                  # Cu, 8 grains, 1e16 m^-2
  python3 make_polycrystal_loops.py --grains 20 --box 300 --density 3e16 --seed 4
  python3 make_polycrystal_loops.py --kinds glide:1,frank:1 --size 25 --type interstitial
  python3 make_polycrystal_loops.py --lattice BCC --element Ta --a0 3.3024 --kinds glide:2,prismatic:1
"""
import argparse
import itertools
import json
import math
import os
import sys

import numpy as np

import exec_input as E

DEFAULTS = {"FCC": dict(element="Cu", a0=3.615, kinds="glide:4,frank:3,prismatic:2,shockley:1"),
            "BCC": dict(element="Ta", a0=3.3024, kinds="glide:3,prismatic:2")}
KINDS = {"FCC": ("glide", "shockley", "frank", "prismatic"), "BCC": ("glide", "prismatic")}

D110 = [np.array(v) for v in sorted(set(itertools.permutations((1, 1, 0))) | set(itertools.permutations((1, -1, 0))))
        if next(x for x in v if x != 0) > 0]
D111 = [np.array(v) for v in itertools.product((1, -1), repeat=3) if v[0] == 1]
D112 = [np.array(v) for v in sorted(set(itertools.permutations((1, 1, 2))) | set(itertools.permutations((1, 1, -2)))
                                    | set(itertools.permutations((1, -1, 2))) | set(itertools.permutations((-1, 1, 2))))]


def reduced(v):
    return np.array(E.reduce_miller([int(round(x)) for x in v]))


# ---------------------------------------------------------------------------- rotations
def random_rotation(rng):
    """Uniform random rotation (Shoemake), as R: rows = crystal directions of the lab axes."""
    u1, u2, u3 = rng.random(3)
    a, b = math.sqrt(1 - u1), math.sqrt(u1)
    x, y, z, w = a * math.sin(2 * math.pi * u2), a * math.cos(2 * math.pi * u2), \
        b * math.sin(2 * math.pi * u3), b * math.cos(2 * math.pi * u3)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def euler_from_R(R):
    """Bunge (phi1, Phi, phi2) in degrees with g = R^T (as DISCO's orientationFromEuler)."""
    g = R.T
    Phi = math.acos(max(-1., min(1., g[2, 2])))
    if abs(math.sin(Phi)) > 1e-9:
        phi1, phi2 = math.atan2(g[2, 0], -g[2, 1]), math.atan2(g[0, 2], g[1, 2])
    else:
        phi1, phi2 = math.atan2(g[0, 1], g[0, 0]), 0.
    return [math.degrees(phi1) % 360, math.degrees(Phi), math.degrees(phi2) % 360]


def R_from_euler(e):
    p1, P, p2 = (math.radians(x) for x in e)
    c1, s1, c, s, c2, s2 = math.cos(p1), math.sin(p1), math.cos(P), math.sin(P), math.cos(p2), math.sin(p2)
    g = np.array([[c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s],
                  [-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s],
                  [s1 * s, -c1 * s, c]])
    return g.T


# ---------------------------------------------------------------------------- crystallography
def fcc_translation(m):
    """b = a0/2 m is an FCC lattice translation (m integer, sum even)."""
    h = np.round(m)
    return np.allclose(m, h, atol=1e-9) and int(abs(h.sum())) % 2 == 0


def shockley_stable(m, n):
    """lattice.hpp shockleyFaultIsStable: d nhat - b is a lattice translation."""
    k = 2 / (math.sqrt(3) * np.linalg.norm(n))
    return fcc_translation(k * np.asarray(n, float) - np.asarray(m, float))


def draw_loop(lat, kind, typ, size, rng):
    """Crystal-frame loop: shape dict + m, n (b = a0/2 m)."""
    sgn = rng.choice((1, -1))
    if lat == "FCC":
        n = D111[rng.integers(4)] * sgn
        if kind == "glide":
            m = [v for v in D110 if v @ n == 0][rng.integers(3)] * rng.choice((1, -1))
            return dict(shape="circular", radius=size), m.astype(float), n
        if kind == "shockley":
            cands = [v / 3. for v in D112 if v @ n == 0]
            cands = [c for c in cands + [-c for c in cands] if shockley_stable(c, n)]
            return dict(shape="circular", radius=size), cands[rng.integers(len(cands))], n
        if kind == "frank":                      # default vertex axis: <110> vertices and sides
            m = (-2. / 3. if typ == "interstitial" else 2. / 3.) * n
            return dict(shape="hexagonal", circumradius=size), m, n
        b = D110[rng.integers(6)] * sgn          # prismatic rhombus: reach = side sqrt(2/3)
        p1, p2 = [p for p in D111 if p @ b == 0]
        d1, d2 = reduced(np.cross(b, p1)), reduced(np.cross(b, p2))
        if d1 @ d2 < 0:
            d2 = -d2
        side = size / math.sqrt(2. / 3.)
        m = (-1. if typ == "interstitial" else 1.) * b
        return dict(shape="parallelogram", half_1=side / 2, half_2=side / 2,
                    dir_1=[int(x) for x in d1], dir_2=[int(x) for x in d2]), m, b
    # BCC
    if kind == "glide":
        m = D111[rng.integers(4)] * sgn
        n = [v for v in D110 if v @ m == 0][rng.integers(3)]
        return dict(shape="circular", radius=size), m.astype(float), n
    b = D111[rng.integers(4)] * sgn              # prismatic hexagon, <112> vertices -> <112> sides
    vdir = next(v for v in D112 if v @ b == 0)
    m = (-1. if typ == "interstitial" else 1.) * b
    return dict(shape="hexagonal", circumradius=size, vertex_dir=[int(x) for x in vdir]), m, b


def perimeter(shape):
    if shape["shape"] == "circular":
        return 2 * math.pi * shape["radius"]
    if shape["shape"] == "hexagonal":
        return 6 * shape["circumradius"]
    return 4 * 2 * shape["half_1"]


def reach(shape):
    if shape["shape"] == "circular":
        return shape["radius"]
    if shape["shape"] == "hexagonal":
        return shape["circumradius"]
    return 2 * shape["half_1"] * math.sqrt(2. / 3.)


def min_image(r, box):
    return r - box * np.round(r / box)


# ---------------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--output", default="polycrystal_loops.yaml")
    p.add_argument("--lattice", default="FCC", choices=["FCC", "BCC"])
    p.add_argument("--element", default=None)
    p.add_argument("--mass", type=float, default=None)
    p.add_argument("--a0", type=float, default=None)
    p.add_argument("--box", type=float, nargs="+", default=[150.], help="Lx [Ly Lz] (A)")
    p.add_argument("--grains", type=int, default=8)
    p.add_argument("--density", type=float, default=1e16, help="nominal line length / volume (m^-2)")
    p.add_argument("--loops", type=int, default=None, help="number of loops instead of --density")
    p.add_argument("--kinds", default=None, help="name:weight list (default: lattice dependent)")
    p.add_argument("--size", type=float, default=20., help="median loop radius / reach (A)")
    p.add_argument("--size-sigma", type=float, default=0.3, help="log-sd of the size (0: all equal)")
    p.add_argument("--type", default="mixed", choices=["interstitial", "vacancy", "mixed"])
    p.add_argument("--gap", type=float, default=10.0, help="minimum distance between loops (A)")
    p.add_argument("--gb-min-distance", type=float, default=None)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--format", default="lammps", choices=["xyz", "lammps", "both", "none"])
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    lat = a.lattice
    a0 = a.a0 or DEFAULTS[lat]["a0"]
    element = a.element or DEFAULTS[lat]["element"]
    box = np.array((a.box * 3)[:3] if len(a.box) == 1 else a.box, float)
    if len(box) != 3:
        sys.exit("--box takes 1 or 3 lengths")

    kinds, weights = [], []
    for item in (a.kinds or DEFAULTS[lat]["kinds"]).split(","):
        name, _, w = item.partition(":")
        if name not in KINDS[lat]:
            sys.exit(f"--kinds: '{name}' is not a {lat} loop kind ({', '.join(KINDS[lat])})")
        kinds.append(name)
        weights.append(float(w or 1))
    weights = np.array(weights) / sum(weights)

    # grains
    grains = []
    for _ in range(a.grains):
        R = random_rotation(rng)
        e = euler_from_R(R)
        grains.append(dict(center=list((rng.random(3) - 0.5) * box), euler=e, R=R_from_euler(e)))
    seeds = np.array([g["center"] for g in grains])

    def grain_of(x):
        return int(np.argmin(np.linalg.norm(min_image(x - seeds, box), axis=1)))

    def cell_margin(x, k):
        """Distance from x to the nearest face of grain k's cell (> 0 inside)."""
        r = min_image(x - seeds[k], box)
        best = np.inf
        for j in range(len(seeds)):
            for sh in itertools.product((-1, 0, 1), repeat=3):
                if j == k and sh == (0, 0, 0):
                    continue
                d = min_image(seeds[j] - seeds[k], box) + np.array(sh) * box
                best = min(best, (d @ d / 2 - r @ d) / np.linalg.norm(d))
        return best

    # loops
    V = box.prod()
    target = a.density * V * 1e-20 if a.loops is None else None
    loops, placed, total = [], [], 0.
    for attempt in range(200000):
        if (a.loops is not None and len(loops) >= a.loops) or (target is not None and total >= target):
            break
        kind = kinds[rng.choice(len(kinds), p=weights)]
        typ = a.type if a.type != "mixed" else ("interstitial", "vacancy")[rng.integers(2)]
        size = a.size * (math.exp(a.size_sigma * rng.standard_normal()) if a.size_sigma > 0 else 1.)
        shape, m, n = draw_loop(lat, kind, typ, size, rng)
        rr = reach(shape)
        if 2 * rr + a.gap > box.min():
            continue
        c = (rng.random(3) - 0.5) * box
        if any(np.linalg.norm(min_image(c - c2, box)) < rr + r2 + a.gap for c2, r2 in placed):
            continue
        k = grain_of(c)
        placed.append((c, rr))
        total += perimeter(shape)
        loops.append(dict(**shape, center=list(c), m=[float(x) for x in m], n=[int(x) for x in n],
                          grain=k, kind=kind, type=typ if kind in ("frank", "prismatic") else None,
                          may_clip=bool(cell_margin(c, k) < rr)))
    else:
        sys.exit("could not reach the requested density: lower --density / --size / --gap")
    if not loops:
        sys.exit("no loop placed")

    d_min = min(E.mid_spacing(l["n"], a0, lat) for l in loops)
    general = E.general_block(a0, (0, 0, 0), 0.33 * d_min, a.format, lattice=lat, element=element, mass=a.mass)
    general["polycrystal"] = dict(box=list(box), gb_min_distance=a.gb_min_distance,
                                  grains=[dict(center=g["center"], euler=g["euler"]) for g in grains])
    counts = {k: sum(l["kind"] == k for l in loops) for k in kinds}
    header = [f"{lat} {element} polycrystal: {len(grains)} grains, {len(loops)} loops "
              f"({', '.join(f'{v} {k}' for k, v in counts.items())}); generated by make_polycrystal_loops.py",
              f"nominal density {total / V * 1e20:.3e} m^-2 (before clipping at the GBs); median size {a.size} A",
              "m, n: crystal frame of the loop's grain; Center: box frame (origin at the box centre)"]
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    E.write_yaml(a.output, general, loops, header)

    side = dict(generator="make_polycrystal_loops.py",
                command=" ".join([os.path.basename(sys.argv[0])] + sys.argv[1:]),
                lattice=lat, element=element, a0=a0, box=list(box),
                nominal_density_m2=total / V * 1e20,
                grains=[dict(index=i + 1, center=g["center"], euler_deg=g["euler"], R=g["R"].tolist())
                        for i, g in enumerate(grains)],
                loops=[dict(index=i + 1, grain=l["grain"] + 1, kind=l["kind"], type=l["type"], shape=l["shape"],
                            size=reach(l), perimeter=perimeter(l), center=l["center"], m=l["m"], n=l["n"],
                            burgers_crystal_a0=[x / 2 for x in l["m"]],
                            burgers_lab_A=list(grains[l["grain"]]["R"] @ (a0 / 2 * np.array(l["m"]))),
                            may_clip=l["may_clip"]) for i, l in enumerate(loops)],
                note="DISCO clips loops at the GBs and skips loops that miss their grain; "
                     "see Poly_..._loops.txt for the final list and density")
    with open(os.path.splitext(a.output)[0] + ".json", "w") as f:
        json.dump(E._jsonable(side), f, indent=1)
    print(f"wrote {a.output}: {len(grains)} grains, {len(loops)} loops "
          f"({', '.join(f'{v} {k}' for k, v in counts.items())}), nominal density {total / V * 1e20:.2e} m^-2, "
          f"{sum(l['may_clip'] for l in loops)} may be clipped")


if __name__ == "__main__":
    main()
