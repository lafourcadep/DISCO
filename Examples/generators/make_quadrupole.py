#!/usr/bin/env python3
"""Write an exec_cuda input for a quadrupole of four glissile loops on two {110} planes,
or for many quadrupoles in one sample (--tiles).

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

Several quadrupoles in one sample (--tiles NX NY NZ): one quadrupole per tile of a regular
grid of cubic tiles, all written to one input file with one JSON sidecar. Each tile takes a
(plane pair, mode) combination from a pool: the --planes pair, or all 15 pairs, or the 3
perpendicular / 12 oblique ones (--tile-pairs); junction, annihilation or both (--mode mixed).
--order cycle takes the combinations in turn with the preferred Burgers assignment;
--order random draws the combination and a valid assignment per tile (--seed). --noise
shifts each tile rigidly by a random vector rounded to whole a0: a lattice translation, so
every loop plane stays on its mid-plane. Tiles have an even edge (--tile-cells, by default
the smallest one holding every tile's quadrupole with a 12 a0 margin, plus room for the
noise), so the tile centres are lattice vectors too. Every loop must stay 2 a0 inside its
tile. Square/rectangular tubes need one Loopon_size for the whole file, so all their tiles
must have the same plane angle (--tile-pairs perpendicular or oblique).

Usage:
  python3 make_quadrupole.py                                   # junction, (110)/(-110), circles
  python3 make_quadrupole.py --mode annihilation --shape square
  python3 make_quadrupole.py --planes 1 1 0 0 1 1 --shape rectangular        # 60 deg, <111> line
  python3 make_quadrupole.py --shape square --d 40
  python3 make_quadrupole.py --shape rectangular --half-m 60     # 60 A half-length along L
  python3 make_quadrupole.py --all-pairs --shape square --outdir quads/
  python3 make_quadrupole.py --keep-offset --cells 40 --units a0 --radius 11 -o quad_junction.yaml
                                  # centres at exactly --d, snapped by exec_cuda
  python3 make_quadrupole.py --tiles 2 2 1                     # 4 identical junction quadrupoles
  python3 make_quadrupole.py --tiles 3 3 3 --tile-pairs all --mode mixed --order random --seed 1 --noise 5
  python3 make_quadrupole.py --tiles 2 2 2 --shape square --tile-pairs perpendicular --mode mixed
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


def design_quadrupole(a, n1, n2, mode, unit_len, rng=None):
    """Loops of one quadrupole centred on the origin (no box yet). The Burgers assignment is
    the preferred one, or a random valid one when rng is given."""
    d110 = E.mid_spacing(n1, a.a0)
    s_req = a.d * unit_len
    # {110} atomic planes pass through the box centre for any Nx, so the mid-plane rounding
    # of the plane offset does not depend on the box size
    s = s_req if a.keep_offset else E.mid_plane_offset(s_req, n1, a.a0, (2, 2, 2))
    L, faces, centres, corners = design_geometry(n1, n2, s)
    combos = search(n1, n2, mode, L, faces, corners)
    if not combos:
        sys.exit(f"no {mode} Burgers assignment exists for planes {list(n1)} / {list(n2)}")
    if rng is None:
        pref = FIGURE4.get((tuple(int(x) for x in n1), tuple(int(x) for x in n2)), {}).get(mode)
        combo = tuple(pref) if pref and tuple(pref) in combos else combos[0]
    else:
        combo = combos[int(rng.integers(len(combos)))]

    Lint = np.cross(n1, n2)
    Lint = Lint // np.gcd.reduce(np.abs(Lint))
    shared = a.shape != "circular"
    h = w = None
    if shared:
        loops, da, h, w = shared_side_loops(a, faces, centres, corners, combo, Lint, unit_len)
    else:
        dims = dict(radius=a.radius * unit_len)
        if not a.size_given:
            dims = auto_scale(faces, centres, combo, dims, L, a.overlap)
        loops = make_loops(faces, centres, combo, dims)
        da = a.loopon_size * unit_len if a.loopon_size else 0.4 * a.a0
    angle = float(np.degrees(np.arccos(abs(unit(n1) @ unit(n2)))))
    return dict(mode=mode, n1=n1, n2=n2, s=s, d110=d110, L=L, Lint=Lint, corners=corners,
                combos=combos, combo=combo, shared=shared, loops=loops, da=da, h=h, w=w, angle=angle)


def list_assignments(Q):
    print(f"{len(Q['combos'])} valid {Q['mode']} assignments for {list(Q['n1'])} / {list(Q['n2'])} (loops 1..4):")
    for cmb in Q["combos"]:
        print("   " + "  ".join(str(list(m)) for m in cmb) + ("   <- used" if cmb == Q["combo"] else ""))


def corner_reactions(Q, a, bcc, cells):
    """Checks the crossing (circles) or the shared sides (tubes) on the centres exec_cuda
    will use in a cubic box of `cells`, and returns the four corner reactions."""
    loops, L, da = Q["loops"], Q["L"], Q["da"]
    snapped = [dict(l, center=E.snap_center(l["center"], l["n"], a.a0, (cells,) * 3)) for l in loops]
    if Q["shared"]:
        check_shared(snapped, L, da)
        h = Q["h"]
        chords = {c: dict(chord_A=[-h, h], chord_B=[-h, h], overlap=2 * h, shared_side=True)
                  for c in CORNERS}
    else:
        chords = check_crossing(snapped, L, a.warn_overlap)
    reactions = []
    for (A, B), c in Q["corners"].items():
        mj = c["signs"][0] * np.array(Q["combo"][A]) + c["signs"][1] * np.array(Q["combo"][B])
        bj = [int(x) // 2 for x in mj]                         # m_J = 2 * (junction b / a0)
        frank = abs(mj @ mj) / 4 < 1.5                         # |b_J|^2 < |b_A|^2 + |b_B|^2 (units a0^2)
        reactions.append(dict(loops=[A + 1, B + 1], line=[int(x) for x in Q["Lint"]], product_a0=bj,
                              kind=product_kind(mj), frank_favourable=bool(frank),
                              expected_junction_type=E.junction_type(bcc, bj) if frank else None,
                              **chords[(A, B)]))
    return reactions


def fmt(v):
    return "[" + " ".join(str(int(x)) for x in v) + "]"


def note_for(shared):
    return ("the four loops share their sides along L: each corner leg carries the product b directly "
            "in the unrelaxed sample" if shared else
            "corner reactions are expected after relaxation; the unrelaxed sample has 4 crossing loops")


def one_quadrupole(a, bcc, n1, n2, unit_len, outpath):
    Q = design_quadrupole(a, n1, n2, a.mode, unit_len)
    if a.list:
        list_assignments(Q)
    loops, shared, Lint, s, angle = Q["loops"], Q["shared"], Q["Lint"], Q["s"], Q["angle"]
    cells = a.cells or E.auto_cells(loops, a.a0)
    E.check_in_box(loops, a.a0, (cells,) * 3)
    reactions = corner_reactions(Q, a, bcc, cells)
    general = E.general_block(a.a0, cells, Q["da"], a.format)

    header = [f"Quadrupole ({a.mode}) of {a.shape} loops on {fmt(n1)} / {fmt(n2)}, generated by make_quadrupole.py",
              f"Prism along L = {fmt(Lint)} ({angle:.0f} deg between the planes); loop planes at "
              f"n.x = +/-{s:.4f} A = {s / Q['d110']:.2f} d(110)",
              "Loop order: 1 = face -n1, 2 = face +n2, 3 = face -n2, 4 = face +n1"]
    if shared:
        header += [f"Closed tube: sides along L (Side) shared exactly by adjacent loops; face half-width "
                   f"{Q['w']:.4f} A = {(Q['w'] / Q['da']) - 0.5:.0f}.5 Loopon_size, half-length along L {Q['h']:.4f} A"]
    header += ["Corner products (both lines along +L), b in units of a0:"]
    header += [f"  {r['loops'][0]}-{r['loops'][1]}: {fmt(r['product_a0'])}  {r['kind']}"
               + ("  (Frank-favourable junction)" if r["frank_favourable"] else "  (no junction)")
               for r in reactions]
    E.write_yaml(outpath, general, loops, header)
    expected = dict(loops=4, corner_reactions=reactions, note=note_for(shared))
    E.write_sidecar(os.path.splitext(outpath)[0] + ".json", bcc, general, loops, expected,
                    extra=dict(quadrupole=dict(mode=a.mode, n1=[int(x) for x in n1], n2=[int(x) for x in n2],
                                               L=[int(x) for x in Lint], angle_deg=angle, plane_offset=s)))
    print(f"{outpath}: {a.mode}, {a.shape}, {fmt(n1)}/{fmt(n2)} (L {fmt(Lint)}, {angle:.0f} deg), "
          f"box {cells}^3, corners " + ", ".join(f"{r['loops'][0]}-{r['loops'][1]} {fmt(r['product_a0'])}"
                                                   for r in reactions))


# =============================================================================
# Several quadrupoles in one sample (--tiles)
# =============================================================================
def tile_pool(a, planes):
    """(n1, n2, mode) combinations the tiles are drawn from."""
    if a.tile_pairs == "planes":
        pairs = [planes]
    else:
        pairs = list(itertools.combinations(planes_110(), 2))
        perp = lambda p: abs(p[0] @ p[1]) == 0
        if a.tile_pairs == "perpendicular":
            pairs = [p for p in pairs if perp(p)]
        elif a.tile_pairs == "oblique":
            pairs = [p for p in pairs if not perp(p)]
    modes = ["junction", "annihilation"] if a.mode == "mixed" else [a.mode]
    return [(n1, n2, m) for n1, n2 in pairs for m in modes]


def tile_shift(a, rng, amp):
    """Rigid random shift of one tile, each component a whole multiple of a0 (a lattice
    vector), so every loop plane stays on its mid-plane. Gaussian draws are cut at 2 sigma."""
    shift = np.zeros(3)
    if amp <= 0:
        return shift
    for ax in a.noise_axes:
        i = "xyz".index(ax)
        x = rng.uniform(-amp, amp) if a.noise_type == "uniform" else float(np.clip(rng.normal(0, amp), -2 * amp, 2 * amp))
        shift[i] = round(x / a.a0) * a.a0
    return shift


def tiled_quadrupoles(a, bcc, planes, unit_len, outpath):
    rng = np.random.default_rng(a.seed)
    NX, NY, NZ = a.tiles
    nt = NX * NY * NZ
    pool = tile_pool(a, planes)
    if a.order == "cycle":
        picks = [pool[t % len(pool)] for t in range(nt)]
    else:
        picks = [pool[int(i)] for i in rng.integers(len(pool), size=nt)]
    Qs = [design_quadrupole(a, n1, n2, mode, unit_len, rng if a.order == "random" else None)
          for n1, n2, mode in picks]

    das = sorted({round(float(Q["da"]), 10) for Q in Qs})
    if len(das) > 1:
        sys.exit(f"square/rectangular tiles need one Loopon_size for the whole sample, but these plane pairs "
                 f"give {das} A: use --tile-pairs perpendicular or oblique (all pairs at the same angle), "
                 f"or circular loops")
    da = Qs[0]["da"]

    amp = a.noise * unit_len
    reach = 2 * amp if a.noise_type == "gaussian" else amp            # largest possible shift
    if a.tile_cells:
        T = a.tile_cells
        if T % 2:
            sys.exit("--tile-cells must be even: tile centres must be lattice vectors so the loop planes "
                     "stay on mid-planes")
    else:
        T = max(E.auto_cells(Q["loops"], a.a0) for Q in Qs) + 2 * int(np.ceil(reach / a.a0))
        T += T % 2
    margin = 2 * a.a0                                                 # loops of two tiles >= 4 a0 apart

    loops, tiles = [], []
    for t, (Q, (i, j, k)) in enumerate(zip(Qs, itertools.product(range(NX), range(NY), range(NZ)))):
        reactions = corner_reactions(Q, a, bcc, T)                    # checked in the tile's own box
        centre = np.array([i - (NX - 1) / 2, j - (NY - 1) / 2, k - (NZ - 1) / 2]) * T * a.a0
        shift = tile_shift(a, rng, amp)
        for l in Q["loops"]:
            ext = np.abs(np.asarray(l["center"], float) + shift) + E.loop_half_extent(l)
            if np.any(ext > 0.5 * T * a.a0 - margin):
                sys.exit(f"tile {t + 1}: a loop reaches {ext.max():.1f} A from the tile centre, more than the "
                         f"half-tile {0.5 * T * a.a0:.1f} A minus a {margin:.1f} A margin: "
                         f"increase --tile-cells or reduce --noise")
            loops.append(dict(l, center=centre + shift + np.asarray(l["center"], float)))
        off = 4 * t
        for r in reactions:
            r["loops"] = [x + off for x in r["loops"]]
        tiles.append(dict(tile=t + 1, grid_index=[i, j, k], centre=centre + shift, shift=shift, mode=Q["mode"],
                          n1=[int(x) for x in Q["n1"]], n2=[int(x) for x in Q["n2"]],
                          L=[int(x) for x in Q["Lint"]], angle_deg=Q["angle"], plane_offset=Q["s"],
                          loops=[off + 1, off + 2, off + 3, off + 4],
                          m=[list(m) for m in Q["combo"]], corner_reactions=reactions))

    cells = (T * NX, T * NY, T * NZ)
    if not a.keep_offset:
        for n, l in enumerate(loops, 1):                              # tiles moved by lattice vectors only
            moved = np.linalg.norm(E.snap_center(l["center"], l["n"], a.a0, cells) - l["center"])
            if moved > 1e-6:
                sys.exit(f"internal error: loop {n} is not on a mid-plane after tiling (moved {moved:.2e} A)")
    general = E.general_block(a.a0, cells, da, a.format)

    nj = sum(r["frank_favourable"] for tl in tiles for r in tl["corner_reactions"])
    header = [f"{nt} quadrupoles ({a.mode}) of {a.shape} loops on a {NX} x {NY} x {NZ} grid, "
              f"generated by make_quadrupole.py",
              f"Tiles of {T}^3 cells; plane pairs: {a.tile_pairs}; order: {a.order}"
              + (f"; seed {a.seed}" if a.seed is not None else ""),
              (f"Rigid tile shifts: {a.noise_type}, amplitude {amp:.4f} A on {a.noise_axes}, rounded to whole a0"
               if amp > 0 else "No tile shifts"),
              f"{4 * nt} loops, {nj} Frank-favourable corners (a<100> junctions), {4 * nt - nj} without junction",
              "Per tile: loops, planes n1 / n2, prism axis L, Burgers m of loops 1..4"]
    header += [f"  tile {tl['tile']:>3}: loops {tl['loops'][0]}-{tl['loops'][3]}, {tl['mode']:<12} "
               f"{fmt(tl['n1'])}/{fmt(tl['n2'])} L {fmt(tl['L'])}  m " + " ".join(fmt(m) for m in tl["m"])
               for tl in tiles]
    E.write_yaml(outpath, general, loops, header)
    expected = dict(loops=4 * nt, quadrupoles=nt, frank_favourable_corners=nj, tiles=tiles,
                    note=note_for(a.shape != "circular"))
    E.write_sidecar(os.path.splitext(outpath)[0] + ".json", bcc, general, loops, expected,
                    extra=dict(tiling=dict(tiles=[NX, NY, NZ], tile_cells=T, tile_pairs=a.tile_pairs,
                                           order=a.order, seed=a.seed, noise=amp, noise_type=a.noise_type,
                                           noise_axes=a.noise_axes, mode=a.mode, shape=a.shape)))
    print(f"{outpath}: {nt} {a.mode} quadrupoles ({a.shape}), {len({(fmt(q['n1']), fmt(q['n2'])) for q in Qs})} "
          f"plane pair(s), tiles {T}^3, box {cells[0]}x{cells[1]}x{cells[2]} "
          f"({2 * cells[0] * cells[1] * cells[2]:,} atoms), {4 * nt} loops, {nj} junction corners")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", default="junction", choices=["junction", "annihilation", "mixed"],
                   help="mixed (with --tiles): each tile is a junction or an annihilation quadrupole")
    p.add_argument("--planes", type=int, nargs=6, metavar="H", default=[1, 1, 0, -1, 1, 0],
                   help="n1 then n2, two different {110} planes (default: 1 1 0 -1 1 0)")
    p.add_argument("--all-pairs", action="store_true", help="write the 15 {110} plane pairs, one file each")
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
    p.add_argument("-o", "--output", help="output .yaml (single pair, or one tiled sample)")
    p.add_argument("--outdir", default=".")
    t = p.add_argument_group("several quadrupoles in one sample")
    t.add_argument("--tiles", type=int, nargs=3, metavar="N",
                   help="NX NY NZ: one quadrupole per tile of a NX x NY x NZ grid, all in one file")
    t.add_argument("--tile-cells", type=int, help="tile edge in cells, even (default: the smallest even cube "
                   "holding every tile's quadrupole with a 12 a0 margin, enlarged for --noise)")
    t.add_argument("--tile-pairs", default="planes", choices=["planes", "all", "perpendicular", "oblique"],
                   help="plane pairs the tiles use: the --planes pair, all 15 {110} pairs, the 3 "
                        "perpendicular ones (L = <100>) or the 12 at 60 deg (L = <111>)")
    t.add_argument("--order", default="cycle", choices=["cycle", "random"],
                   help="cycle: tiles take the (pair, mode) combinations in turn with the preferred Burgers "
                        "assignment; random: random combination and random valid assignment per tile")
    t.add_argument("--noise", type=float, default=0.0,
                   help="amplitude of a rigid random shift of each tile (std for gaussian, cut at 2 std; "
                        "half-width for uniform); shifts are rounded to whole a0")
    t.add_argument("--noise-type", default="gaussian", choices=["gaussian", "uniform"])
    t.add_argument("--noise-axes", default="xyz", help="axes the shifts act on, any subset of xyz")
    t.add_argument("--seed", type=int, help="random seed (--order random, --noise)")
    a = p.parse_args()

    unit_len = a.a0 if a.units == "a0" else 1.0
    if a.shape == "circular" and a.half_m is not None:
        sys.exit("--half-m is for rectangular loops")
    if a.shape != "circular" and a.radius is not None:
        sys.exit("--radius is for circular loops; square/rectangular loops span the face set by --d")
    if a.shape != "circular" and a.keep_offset:
        sys.exit("--keep-offset is for circular loops: shared sides need the faces on mid-planes")
    if a.tiles:
        if min(a.tiles) < 1:
            sys.exit("--tiles: every count must be >= 1")
        for opt, bad in (("--all-pairs", a.all_pairs), ("--cells", a.cells), ("--list", a.list),
                         ("--keep-offset", a.keep_offset)):
            if bad:
                sys.exit(f"{opt} is for a single quadrupole; with --tiles use --tile-cells / --tile-pairs")
        if set(a.noise_axes) - set("xyz") or not a.noise_axes:
            sys.exit("--noise-axes: any subset of xyz")
        if a.noise < 0:
            sys.exit("--noise must be >= 0")
    else:
        if a.mode == "mixed":
            sys.exit("--mode mixed needs --tiles")
        tiled_opts = [o for o, v in (("--tile-cells", a.tile_cells), ("--noise", a.noise),
                                     ("--seed", a.seed is not None), ("--order random", a.order == "random"),
                                     ("--tile-pairs", a.tile_pairs != "planes")) if v]
        if tiled_opts:
            sys.exit(f"{', '.join(tiled_opts)}: only with --tiles")
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
    if a.tiles:
        n1, n2 = pairs[0]
        what = f"{tag(n1)}_{tag(n2)}" if a.tile_pairs == "planes" else a.tile_pairs
        out = a.output or os.path.join(
            a.outdir, f"quad_tiles_{'x'.join(map(str, a.tiles))}_{a.mode}_{what}_{SHAPE_TAG[a.shape]}.yaml")
        tiled_quadrupoles(a, bcc, (n1, n2), unit_len, out)
        return
    for n1, n2 in pairs:
        out = a.output or os.path.join(a.outdir, f"quad_{a.mode}_{tag(n1)}_{tag(n2)}_{SHAPE_TAG[a.shape]}.yaml")
        one_quadrupole(a, bcc, n1, n2, unit_len, out)


if __name__ == "__main__":
    main()
