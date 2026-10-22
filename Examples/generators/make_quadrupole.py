#!/usr/bin/env python3
"""Write an exec_cuda input for a quadrupole of four glissile loops on two {110} planes.

Junction and annihilation quadrupoles on any pair of non-parallel {110} planes n1, n2,
with circular, square or rectangular loops.

Geometry: the four loops lie on the planes n1.x = -/+s and n2.x = +/-s, the faces of a
prism along L = n1 x n2 (seen along L: a square if n1 is normal to n2, which is 3 pairs
with L = <100>; a 60 deg rhombus otherwise, which is 12 pairs with L = <111>). Each loop is
centred on the middle of its face.
  square / rectangular  the loops are the faces of a closed tube: each spans exactly its face
                width, with two sides along L (exec_cuda key Side = L) that coincide with the
                neighbours' sides on the four corner lines. The shared legs carry b_A +/- b_B
                directly. Exact coincidence needs the loopon tile edge on the corner line, so
                Loopon_size is tuned (as in make_prismatic_dipole.py): face half-width
                w = (k + 1/2) Loopon_size. Square: half-length along L = w; rectangle:
                --half-m along L (default 14/11 w), rounded to (k' + 1/2) Loopon_size.
  circular      adjacent loops cross on the corner lines (they overlap instead of touching;
                exact tangency gives DXA artefacts).
Loop order: 1 = face -n1, 2 = face +n2, 3 = face -n2, 4 = face +n1.
The planes are put on mid-planes between atomic planes (s rounded down to one), so
exec_cuda's mid-plane snap does not move the loops and the prism stays symmetric.

Burgers vectors: each loop takes a 1/2<111> vector b = a0/2 m of its plane (b parallel to
L excluded: it would be pure screw on both planes). With the line sense right-handed about
the written n, the corner product with both lines oriented along +L is
m_J = s_A m_A + s_B m_B, where s = sign((n x e).L) and e points from the loop centre
to the corner line. The assignment is searched:
  junction      all 4 corner products are a<100> (Frank-favourable binary junction)
  annihilation  opposite faces carry +/-b (two dipoles), all 4 corner products are <110>
                (Frank-unfavourable: no junction, the dipoles annihilate by cross-slip)
Frank's rule is also applied, in both cases, at each loop-loop corner.

Loop size: square and rectangular loops are sized by --d (face width) as above. A circle,
without --radius, is scaled to the smallest radius, rounded up to 0.1 A, whose chords on
every corner line overlap by --overlap (20 A); an explicit radius is checked, and it is
rejected if two adjacent loops do not cross. The box is the smallest even cube keeping a
12 a0 margin, unless --cells is given.

Usage:
  python3 make_quadrupole.py                                   # junction, (110)/(-110), circles
  python3 make_quadrupole.py --mode annihilation --shape square
  python3 make_quadrupole.py --planes 1 1 0 0 1 1 --shape rectangular        # 60 deg, <111> line
  python3 make_quadrupole.py --shape square --d 40
  python3 make_quadrupole.py --shape rectangular --half-m 60     # 60 A half-length along L
  python3 make_quadrupole.py --all-pairs --shape square --outdir quads/
  python3 make_quadrupole.py --keep-offset --cells 40 --units a0 --radius 11 -o quad_junction.yaml
                                  # centres at exactly --d, snapped by exec_cuda
"""
import argparse
import itertools
import os
import sys

import numpy as np

import exec_input as E

# hand-built assignments for (110)/(-110), kept as the preferred solution
FIGURE4 = {((1, 1, 0), (-1, 1, 0)): {
    "junction": [(1, -1, -1), (-1, -1, 1), (-1, -1, -1), (1, -1, 1)],
    "annihilation": [(1, -1, -1), (-1, -1, -1), (1, 1, 1), (-1, 1, 1)]}}
CORNERS = [(0, 1), (0, 2), (3, 1), (3, 2)]      # adjacent loops (0-based)
SHAPE_TAG = {"circular": "circ", "square": "sq", "rectangular": "rect"}


def planes_110():
    """The six {110} planes, one sign each, in a fixed order."""
    out = []
    for n in [(1, 1, 0), (1, 0, 1), (0, 1, 1), (-1, 1, 0), (-1, 0, 1), (0, -1, 1)]:
        out.append(np.array(n))
    return out


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def design_geometry(n1, n2, s):
    """Face centres and corner points for plane offset s (A)."""
    L = unit(np.cross(n1, n2))
    nh = [unit(n1), unit(n2)]
    c = nh[0] @ nh[1]
    faces = [(n1, -1, 0), (n2, +1, 1), (n2, -1, 1), (n1, +1, 0)]   # (written n, side, family)
    centres = []
    for _, sg, i in faces:
        a = sg * s / (1 - c * c)
        centres.append(a * nh[i] - c * a * nh[1 - i])
    corners = {}
    for A, B in CORNERS:
        M = np.array([nh[faces[A][2]], nh[faces[B][2]], L])
        K = np.linalg.solve(M, [faces[A][1] * s, faces[B][1] * s, 0.])
        signs = [int(np.sign(np.cross(unit(faces[X][0]), unit(K - centres[X])) @ L)) for X in (A, B)]
        corners[(A, B)] = dict(point=K, signs=signs)
    return L, faces, centres, corners


def product_kind(v):
    a = sorted(abs(int(x)) for x in v)
    return "a<100>" if a == [0, 0, 2] else "a<110>" if a == [0, 2, 2] else "other"


def search(n1, n2, mode, L, faces, corners):
    """All Burgers assignments (tuples of 4 integer m) that satisfy the mode."""
    m111 = [np.array(v) for v in itertools.product((1, -1), repeat=3)]
    cands = [[m for m in m111 if m @ f[0] == 0 and abs(abs(unit(m) @ L) - 1) > 1e-9] for f in faces]
    found = []
    for combo in itertools.product(*cands):
        kinds = [product_kind(c["signs"][0] * combo[A] + c["signs"][1] * combo[B])
                 for (A, B), c in corners.items()]
        if mode == "junction" and all(k == "a<100>" for k in kinds):
            found.append(combo)
        if (mode == "annihilation" and all(k == "a<110>" for k in kinds)
                and np.array_equal(combo[3], -combo[0]) and np.array_equal(combo[2], -combo[1])):
            found.append(combo)
    return [tuple(tuple(int(x) for x in m) for m in combo) for combo in found]


def make_loops(faces, centres, combo, dims):
    """Circular loops (crossing design)."""
    return [dict(shape="circular", center=c, m=list(m), n=[int(x) for x in n], **dims)
            for (n, _, _), c, m in zip(faces, centres, combo)]


def corner_chords(loops, L):
    """For each adjacent pair: the corner line (intersection of the two loop planes as they
    are, after any snap) and the chord each loop cuts on it. Chord = None if the loop misses."""
    out = {}
    for (A, B) in CORNERS:
        nA, nB = unit(loops[A]["n"]), unit(loops[B]["n"])
        K = np.linalg.solve(np.array([nA, nB, L]),
                            [nA @ np.asarray(loops[A]["center"]), nB @ np.asarray(loops[B]["center"]), 0.])
        out[(A, B)] = (K, [E.chord_on_line(loops[X], K, L) for X in (A, B)])
    return out


def overlap_of(chords):
    if None in chords:
        return -np.inf
    return min(chords[0][1], chords[1][1]) - max(chords[0][0], chords[1][0])


def min_overlap(loops, L):
    return min(overlap_of(ch) for _, ch in corner_chords(loops, L).values())


def check_crossing(loops, L, warn_overlap):
    """Each adjacent pair must both cut the corner line, with overlapping chords."""
    info = {}
    for (A, B), (K, chords) in corner_chords(loops, L).items():
        for X, ch in zip((A, B), chords):
            if ch is None:
                c = np.asarray(loops[X]["center"])
                dist, sup = np.linalg.norm(K - c), loops[X]["radius"]
                sys.exit(f"loop {X + 1} does not reach the corner line it shares with loop "
                         f"{(B if X == A else A) + 1}: centre-to-line distance {dist:.2f} A, loop reaches "
                         f"{sup:.2f} A. Enlarge the loops (omit the size options to size them "
                         f"automatically) or reduce --d.")
        ov = overlap_of(chords)
        if ov <= 0:
            sys.exit(f"loops {A + 1} and {B + 1} cut their corner line on disjoint chords, so they do not "
                     f"cross: enlarge the loops (omit the size options to size them automatically)")
        if ov < warn_overlap:
            print(f"WARNING: loops {A + 1}-{B + 1}: chord overlap only {ov:.2f} A on the corner line "
                  f"(near tangency, DXA may leave dangling ends)")
        info[(A, B)] = dict(chord_A=[float(x) for x in chords[0]], chord_B=[float(x) for x in chords[1]],
                            overlap=float(ov))
    return info


def auto_scale(faces, centres, combo, dims, L, target):
    """Smallest scale of the default circle (radius rounded up to 0.1 A) giving
    a chord overlap >= target at every corner."""
    f = lambda lam: min_overlap(make_loops(faces, centres, combo,
                                           {k: v * lam for k, v in dims.items()}), L)
    lo, hi = 0.05, 1.0
    while f(hi) < target:
        hi *= 2
        if hi > 64:
            sys.exit("automatic loop sizing failed")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if f(mid) >= target else (mid, hi)
    return {k: round(float(np.ceil(v * hi / 0.1 - 1e-9)) * 0.1, 6) for k, v in dims.items()}


def shared_side_loops(a, faces, centres, corners, combo, Lint, unit_len):
    """Square/rectangular loops tiling the four prism faces, sides along L shared exactly.
    Returns (loops, Loopon_size, half-length along L, face half-width)."""
    w = {np.linalg.norm(c["point"] - centres[X]) for (A, B), c in corners.items() for X in (A, B)}
    if max(w) - min(w) > 1e-9:
        sys.exit(f"internal error: unequal face widths {sorted(w)}")
    w = max(w)
    target = a.loopon_size * unit_len if a.loopon_size else 0.4 * a.a0
    k_p = max(1, round(w / target - 0.5))
    da = w / (k_p + 0.5)                         # tile edge (k_p + 1/2) da on the corner line
    along = w if a.shape == "square" else (a.half_m * unit_len if a.half_m is not None else 14. / 11. * w)
    k_l = max(1, round(along / da - 0.5))
    # a quarter step past the last grid row: no rounding loss in exec_cuda's |u| <= Half test
    loops = [dict(shape="rectangular", half_m=(k_l + 0.25) * da, half_p=(k_p + 0.25) * da,
                  corner_radius=0., side=[int(x) for x in Lint], center=c, m=list(m),
                  n=[int(x) for x in n]) for (n, _, _), c, m in zip(faces, centres, combo)]
    return loops, da, (k_l + 0.5) * da, w


def check_shared(loops, L, da):
    """After exec_cuda's snap: every corner line lies on both loops' tile edge, (k + 1/2) da
    from the centre across the face, and both loops span the same interval along L."""
    for (A, B), (K, _) in corner_chords(loops, L).items():
        spans = []
        for X in (A, B):
            u, v, _ = E.loop_axes(loops[X])
            r = K - np.asarray(loops[X]["center"], float)
            edge = (np.floor(loops[X]["half_p"] / da) + 0.5) * da
            if abs(r @ u) > 1e-6 or abs(abs(r @ v) - edge) > 1e-6:
                sys.exit(f"loops {A + 1}-{B + 1}: corner line is not on loop {X + 1}'s edge "
                         f"(along L {r @ u:.2e} A, across {abs(r @ v):.6f} A vs edge {edge:.6f} A)")
            spans.append(loops[X]["half_m"])
        if abs(spans[0] - spans[1]) > 1e-9:
            sys.exit(f"loops {A + 1}-{B + 1}: different lengths along L")


def one_quadrupole(a, bcc, n1, n2, unit_len, outpath):
    d110 = E.mid_spacing(n1, a.a0)
    s_req = a.d * unit_len
    # {110} atomic planes pass through the box centre for any Nx, so the mid-plane rounding
    # of the plane offset does not depend on the box size
    s = s_req if a.keep_offset else E.mid_plane_offset(s_req, n1, a.a0, (2, 2, 2))
    L, faces, centres, corners = design_geometry(n1, n2, s)
    combos = search(n1, n2, a.mode, L, faces, corners)
    if not combos:
        sys.exit(f"no {a.mode} Burgers assignment exists for planes {list(n1)} / {list(n2)}")
    pref = FIGURE4.get((tuple(int(x) for x in n1), tuple(int(x) for x in n2)), {}).get(a.mode)
    combo = tuple(pref) if pref and tuple(pref) in combos else combos[0]
    if a.list:
        print(f"{len(combos)} valid {a.mode} assignments for {list(n1)} / {list(n2)} (loops 1..4):")
        for cmb in combos:
            print("   " + "  ".join(str(list(m)) for m in cmb) + ("   <- used" if cmb == combo else ""))

    Lint = np.cross(n1, n2)
    Lint = Lint // np.gcd.reduce(np.abs(Lint))
    shared = a.shape != "circular"
    if shared:
        loops, da, h, w = shared_side_loops(a, faces, centres, corners, combo, Lint, unit_len)
    else:
        dims = dict(radius=a.radius * unit_len)
        if not a.size_given:
            dims = auto_scale(faces, centres, combo, dims, L, a.overlap)
        loops = make_loops(faces, centres, combo, dims)
        da = a.loopon_size * unit_len if a.loopon_size else 0.4 * a.a0
    cells = a.cells or E.auto_cells(loops, a.a0)
    E.check_in_box(loops, a.a0, (cells,) * 3)
    # crossing / sharing is checked on the centres exec_cuda will actually use
    snapped = [dict(l, center=E.snap_center(l["center"], l["n"], a.a0, (cells,) * 3)) for l in loops]
    if shared:
        check_shared(snapped, L, da)
        chords = {c: dict(chord_A=[-h, h], chord_B=[-h, h], overlap=2 * h, shared_side=True)
                  for c in CORNERS}
    else:
        chords = check_crossing(snapped, L, a.warn_overlap)

    general = E.general_block(a.a0, cells, da, a.format)
    angle = float(np.degrees(np.arccos(abs(unit(n1) @ unit(n2)))))
    reactions = []
    for (A, B), c in corners.items():
        mj = c["signs"][0] * np.array(combo[A]) + c["signs"][1] * np.array(combo[B])
        bj = [int(x) // 2 for x in mj]                         # m_J = 2 * (junction b / a0)
        frank = abs(mj @ mj) / 4 < 1.5                         # |b_J|^2 < |b_A|^2 + |b_B|^2 (units a0^2)
        reactions.append(dict(loops=[A + 1, B + 1], line=[int(x) for x in Lint], product_a0=bj,
                              kind=product_kind(mj), frank_favourable=bool(frank),
                              expected_junction_type=E.junction_type(bcc, bj) if frank else None,
                              **chords[(A, B)]))

    fmt = lambda v: "[" + " ".join(str(int(x)) for x in v) + "]"
    header = [f"Quadrupole ({a.mode}) of {a.shape} loops on {fmt(n1)} / {fmt(n2)}, generated by make_quadrupole.py",
              f"Prism along L = {fmt(Lint)} ({angle:.0f} deg between the planes); loop planes at "
              f"n.x = +/-{s:.4f} A = {s / d110:.2f} d(110)",
              "Loop order: 1 = face -n1, 2 = face +n2, 3 = face -n2, 4 = face +n1"]
    if shared:
        header += [f"Closed tube: sides along L (Side) shared exactly by adjacent loops; face half-width "
                   f"{w:.4f} A = {(w / da) - 0.5:.0f}.5 Loopon_size, half-length along L {h:.4f} A"]
    header += ["Corner products (both lines along +L), b in units of a0:"]
    header += [f"  {r['loops'][0]}-{r['loops'][1]}: {fmt(r['product_a0'])}  {r['kind']}"
               + ("  (Frank-favourable junction)" if r["frank_favourable"] else "  (no junction)")
               for r in reactions]
    E.write_yaml(outpath, general, loops, header)
    note = ("the four loops share their sides along L: each corner leg carries the product b directly "
            "in the unrelaxed sample" if shared else
            "corner reactions are expected after relaxation; the unrelaxed sample has 4 crossing loops")
    expected = dict(loops=4, corner_reactions=reactions, note=note)
    E.write_sidecar(os.path.splitext(outpath)[0] + ".json", bcc, general, loops, expected,
                    extra=dict(quadrupole=dict(mode=a.mode, n1=[int(x) for x in n1], n2=[int(x) for x in n2],
                                               L=[int(x) for x in Lint], angle_deg=angle, plane_offset=s)))
    print(f"{outpath}: {a.mode}, {a.shape}, {fmt(n1)}/{fmt(n2)} (L {fmt(Lint)}, {angle:.0f} deg), "
          f"box {cells}^3, corners " + ", ".join(f"{r['loops'][0]}-{r['loops'][1]} {fmt(r['product_a0'])}"
                                                   for r in reactions))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", default="junction", choices=["junction", "annihilation"])
    p.add_argument("--planes", type=int, nargs=6, metavar="H", default=[1, 1, 0, -1, 1, 0],
                   help="n1 then n2, two different {110} planes (default: 1 1 0 -1 1 0)")
    p.add_argument("--all-pairs", action="store_true", help="write the 15 {110} plane pairs")
    p.add_argument("--shape", default="circular", choices=["circular", "square", "rectangular"])
    p.add_argument("--units", default="A", choices=["A", "a0"], help="unit of every length option")
    p.add_argument("--radius", type=float, help="circular: radius (default: automatic, see --overlap)")
    p.add_argument("--half-m", type=float, help="rectangular: half-length along the prism axis L "
                   "(default 14/11 of the face half-width; rounded to the loopon grid)")
    p.add_argument("--d", type=float, help="distance from the prism axis to each loop plane "
                   "(default 7.5*sqrt(2) a0 = 35.03 A, rounded down to a mid-plane)")
    p.add_argument("--keep-offset", action="store_true",
                   help="write the centres at exactly --d and let exec_cuda snap each loop")
    p.add_argument("--loopon-size", type=float, help="loopon grid step (default 0.4 a0 = 1.32096 A); "
                   "square/rectangular: tuned close to it so the shared sides coincide")
    p.add_argument("--cells", type=int, help="Nx = Ny = Nz; default: smallest even cube with a 12 a0 margin")
    p.add_argument("--a0", type=float, default=E.A0)
    p.add_argument("--format", default="xyz", choices=["xyz", "lammps", "both", "none"])
    p.add_argument("--overlap", type=float, default=20.0,
                   help="automatic sizing (no size option given): smallest loops whose chords on each "
                        "corner line overlap by at least this (A)")
    p.add_argument("--warn-overlap", type=float, default=5.0,
                   help="warn when two loops' chords on their corner line overlap by less than this (A)")
    p.add_argument("--list", action="store_true", help="print every valid Burgers assignment")
    p.add_argument("--bcc", default=E.DEFAULT_BCC)
    p.add_argument("-o", "--output", help="output .yaml (single pair only)")
    p.add_argument("--outdir", default=".")
    a = p.parse_args()

    unit_len = a.a0 if a.units == "a0" else 1.0
    if a.shape == "circular" and a.half_m is not None:
        sys.exit("--half-m is for rectangular loops")
    if a.shape != "circular" and a.radius is not None:
        sys.exit("--radius is for circular loops; square/rectangular loops span the face set by --d")
    if a.shape != "circular" and a.keep_offset:
        sys.exit("--keep-offset is for circular loops: shared sides need the faces on mid-planes")
    a.size_given = a.radius is not None
    dflt = lambda x_a0: x_a0 * a.a0 / unit_len
    a.radius = a.radius if a.radius is not None else dflt(11.)
    a.d = a.d if a.d is not None else dflt(7.5 * np.sqrt(2))

    bcc = E.load_bcc(a.bcc)
    if a.all_pairs:
        if a.output:
            sys.exit("-o is for a single pair; use --outdir with --all-pairs")
        pairs = list(itertools.combinations(planes_110(), 2))
    else:
        n1, n2 = np.array(a.planes[:3]), np.array(a.planes[3:])
        for n in (n1, n2):
            if sorted(abs(int(x)) for x in n) != [0, 1, 1]:
                sys.exit(f"{list(n)} is not a {{110}} plane")
        if np.linalg.norm(np.cross(n1, n2)) == 0:
            sys.exit("n1 and n2 are parallel")
        pairs = [(n1, n2)]
    os.makedirs(a.outdir, exist_ok=True)
    tag = lambda n: "".join(str(int(x)) for x in n)
    for n1, n2 in pairs:
        out = a.output or os.path.join(a.outdir, f"quad_{a.mode}_{tag(n1)}_{tag(n2)}_{SHAPE_TAG[a.shape]}.yaml")
        one_quadrupole(a, bcc, n1, n2, unit_len, out)


if __name__ == "__main__":
    main()
