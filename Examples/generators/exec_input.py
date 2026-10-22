"""Shared helpers for the exec_cuda input generators in this folder.

- BCC slip systems read from Parameters/BCC.txt, in DANCE's own order (index 1..48),
  so the slip-system index written in the ground-truth sidecars is the one DANCE reports.
- Python copies of the exec_cuda geometry rules (loop.cpp): mid-plane snap of the loop
  centre, in-plane axes m, p = m x n, and the loopon grid (to predict loopon counts).
- Writers for the YAML input (exact key set accepted by input.cpp) and the JSON sidecar.

Every length is in Angstrom, like in the exec_cuda input.
"""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BCC = os.path.normpath(os.path.join(HERE, "..", "..", "Parameters", "BCC.txt"))
A0 = 3.3024


# =============================================================================
# Slip systems (Parameters/BCC.txt)
# =============================================================================
def plane_family(n):
    """'110', '112' or '123' from integer Miller indices."""
    key = "".join(str(c) for c in sorted(abs(x) for x in reduce_miller(n)))
    return {"011": "110"}.get(key, key)


def load_bcc(path=DEFAULT_BCC):
    """Parse BCC.txt. Returns dict with
    'systems': list of dicts (index, type, b, n, family), 1-based index in file order;
    'junctions': {type: b} for the junction (a<100>) Burgers types."""
    tokens = open(path).read().split()
    i, section, n_std = 0, "standard", None
    systems, junctions, per_type = [], {}, None
    while i < len(tokens):
        t = tokens[i]; i += 1
        if t == "STANDARD":
            n_std = int(tokens[i]); i += 1
        elif t == "TYPES":
            n_types = int(tokens[i]); i += 1
            if section == "standard":
                per_type = n_std // n_types
        elif t == "JUNCTION":
            section = "junction"
        elif t == "TYPE":
            ty = int(tokens[i]); i += 1
            assert tokens[i] == "BURGERS"; i += 1
            b = [float(x) for x in tokens[i:i + 3]]; i += 3
            if section == "junction":
                junctions[ty] = b
                continue
            assert tokens[i] == "SLIP_SYSTEMS"; i += 1
            for _ in range(per_type):
                n = [int(float(x)) for x in tokens[i:i + 3]]; i += 3
                m = [int(round(2 * c)) for c in b]       # b = a0/2 m
                assert sum(x * y for x, y in zip(m, n)) == 0, f"BCC.txt: b.n != 0 for type {ty} n {n}"
                systems.append(dict(index=len(systems) + 1, type=ty, m=m, n=n, family=plane_family(n)))
    return dict(systems=systems, junctions=junctions)


def sign_test(v):
    """DANCE's utils.cpp sign_test: first component with |x| > 1e-3 made positive."""
    for c in v:
        if abs(c) > 1e-3:
            return [x if c > 0 else -x for x in v]
    return list(v)


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def expected_type_and_sys(bcc, m, n):
    """Ground truth (Burgers type, slip-system index) as DANCE should report it,
    or (type, None) if n is not one of that type's slip planes."""
    bh = _unit(sign_test(list(_unit(m))))
    ty = max({s["type"] for s in bcc["systems"]},
             key=lambda t: abs(bh @ _unit(next(s["m"] for s in bcc["systems"] if s["type"] == t))))
    bref = next(s["m"] for s in bcc["systems"] if s["type"] == ty)
    if abs(abs(bh @ _unit(bref)) - 1) > 1e-6:
        return None, None
    nh = _unit(n)
    for s in bcc["systems"]:
        if s["type"] == ty and abs(abs(nh @ _unit(s["n"])) - 1) < 1e-6:
            return ty, s["index"]
    return ty, None


def junction_type(bcc, m):
    """BCC.txt junction type (5, 6, 7) whose Burgers direction is parallel to m, else None."""
    mh = _unit(m)
    for ty, b in bcc["junctions"].items():
        if abs(abs(mh @ _unit(b)) - 1) < 1e-9:
            return ty
    return None


# =============================================================================
# exec_cuda geometry (loop.cpp)
# =============================================================================
def reduce_miller(n):
    n = [int(x) for x in n]
    g = math.gcd(math.gcd(abs(n[0]), abs(n[1])), abs(n[2]))
    return [x // g for x in n]


def mid_spacing(n, a0=A0):
    """Atomic-plane spacing of the family n (as snapToMidPlane in loop.cpp)."""
    h = reduce_miller(n)
    norm = math.sqrt(sum(x * x for x in h))
    return a0 / norm if sum(h) % 2 == 0 else a0 / (2 * norm)


def snap_center(center, n, a0, cells):
    """Centre moved along n onto the mid-plane, exactly as snapToMidPlane in loop.cpp."""
    h = reduce_miller(n)
    nh = _unit(h)
    sp = mid_spacing(h, a0)
    s0 = -0.5 * a0 * np.asarray(cells, float) @ nh
    s = np.asarray(center, float) @ nh
    s_mid = s0 + (math.floor((s - s0) / sp + 1e-9) + 0.5) * sp
    return np.asarray(center, float) + (s_mid - s) * nh


def mid_plane_offset(s, n, a0, cells):
    """Distance s >= 0 along n-hat moved down to the nearest mid-plane between atomic
    planes, so that faces at +s and -s are both left in place by snapToMidPlane."""
    sp = mid_spacing(n, a0)
    ph = (-0.5 * a0 * np.asarray(cells, float) @ _unit(reduce_miller(n))) % sp
    s_mid = ph + (math.floor((s - ph) / sp - 0.5 + 1e-9) + 0.5) * sp
    if abs(((-s_mid - ph) / sp - 0.5) - round((-s_mid - ph) / sp - 0.5)) > 1e-6:
        sys.exit(f"planes {n}: +s and -s cannot both be mid-planes with cells = {cells}; use even cells")
    return s_mid


def in_plane_axes(m, n, side=None):
    """exec_cuda's in-plane loopon-grid axes (loop.cpp): u-hat, the Half_m axis, and
    v-hat = u-hat x n-hat, the Half_p axis. u is m (screw sides) unless a rectangular
    loop gives Side."""
    uh, nh = _unit(m if side is None else side), _unit(n)
    return uh, np.cross(uh, nh), nh


def loop_axes(loop):
    """in_plane_axes of a loop dict (uses its optional 'side'). For a parallelogram: the
    unit vectors of Dir_1, Dir_2 (oblique grid axes) and n-hat."""
    if loop["shape"] == "parallelogram":
        return _unit(loop["dir_1"]), _unit(loop["dir_2"]), _unit(loop["n"])
    return in_plane_axes(loop["m"], loop["n"], loop.get("side"))


def _halves(loop):
    """Grid half-extents (along u, along v) of a non-circular loop."""
    if loop["shape"] == "parallelogram":
        return loop["half_1"], loop["half_2"]
    return loop["half_m"], loop["half_p"]


def inside_2d(loop, u, v, eps=1e-9):
    """Point (u, v) inside the loop contour, with loop.cpp's boundary rules."""
    if loop["shape"] == "circular":
        return math.hypot(u, v) < loop["radius"] - eps
    if loop["shape"] == "parallelogram":        # u, v: oblique coordinates
        return abs(u) <= loop["half_1"] + eps and abs(v) <= loop["half_2"] + eps
    rc = loop["corner_radius"]
    du = max(abs(u) - max(loop["half_m"] - rc, 0.), 0.)
    dv = max(abs(v) - max(loop["half_p"] - rc, 0.), 0.)
    return du * du + dv * dv <= rc * rc + eps


def count_loopons(loop, da):
    """Number of loopons exec_cuda will put in this loop (generateLoopons)."""
    ext_m, ext_p = (loop["radius"],) * 2 if loop["shape"] == "circular" else _halves(loop)
    km, kp = math.ceil(ext_m / da) + 1, math.ceil(ext_p / da) + 1
    return sum(inside_2d(loop, i * da, j * da)
               for i in range(-km, km + 1) for j in range(-kp, kp + 1))


def loop_half_extent(loop):
    """Per-axis half-extent of the loop around its centre (A)."""
    mh, ph, _ = loop_axes(loop)
    if loop["shape"] == "circular":
        return loop["radius"] * np.hypot(mh, ph)
    if loop["shape"] == "parallelogram":
        return loop["half_1"] * np.abs(mh) + loop["half_2"] * np.abs(ph)
    rc = min(loop["corner_radius"], loop["half_m"], loop["half_p"])
    return (loop["half_m"] - rc) * np.abs(mh) + (loop["half_p"] - rc) * np.abs(ph) + rc * np.hypot(mh, ph)


def auto_cells(loops, a0=A0, margin=12 * A0):
    """Smallest even cube Nx = Ny = Nz keeping `margin` A between every loop and the box face."""
    reach = max(np.max(np.abs(np.asarray(l["center"], float)) + loop_half_extent(l)) for l in loops)
    n = math.ceil(2 * (reach + margin) / a0)
    return n + n % 2


def check_in_box(loops, a0, cells, margin=5.0):
    half = 0.5 * a0 * np.asarray(cells, float)
    for i, l in enumerate(loops, 1):
        reach = np.abs(np.asarray(l["center"], float)) + loop_half_extent(l)
        if np.any(reach > half - margin):
            sys.exit(f"loop {i} reaches {reach.max():.1f} A but the half-box is {half.min():.1f} A: "
                     f"increase --cells")


def _clip_box(q, d, hu, hv):
    """t-interval of q + t d inside |u| <= hu, |v| <= hv (Liang-Barsky), or None."""
    t0, t1 = -np.inf, np.inf
    for qi, di, h in ((q[0], d[0], hu), (q[1], d[1], hv)):
        if abs(di) < 1e-15:
            if abs(qi) > h:
                return None
            continue
        a, b = sorted(((-h - qi) / di, (h - qi) / di))
        t0, t1 = max(t0, a), min(t1, b)
    return (t0, t1) if t0 <= t1 else None


def _clip_disk(q, d, c, r):
    """t-interval of q + t d inside the disk |x - c| <= r, or None (d is a unit vector)."""
    w = q - c
    b, cc = w @ d, w @ w - r * r
    disc = b * b - cc
    if disc < 0:
        return None
    return (-b - math.sqrt(disc), -b + math.sqrt(disc))


def chord_on_line(loop, point, direction):
    """Interval [t0, t1] of the line point + t*direction (in the loop plane) inside the
    loop contour, or None. Exact: the rounded rectangle is the union of two cross-shaped
    boxes and four corner disks, and the union of their chords is the chord."""
    mh, ph, nh = loop_axes(loop)
    q3 = np.asarray(point, float) - np.asarray(loop["center"], float)
    d3 = _unit(direction)
    if abs(q3 @ nh) > 1e-6 or abs(d3 @ nh) > 1e-9:
        raise ValueError("chord_on_line: line not in the loop plane")
    if loop["shape"] == "parallelogram":
        # oblique coordinates: solve x = a u + b v in the plane (t stays the length along d3)
        B = np.array([mh, ph]).T
        q = np.linalg.lstsq(B, q3, rcond=None)[0]
        d = np.linalg.lstsq(B, d3, rcond=None)[0]
        return _clip_box(q, d, loop["half_1"], loop["half_2"])
    q = np.array([q3 @ mh, q3 @ ph])
    d = np.array([d3 @ mh, d3 @ ph])
    d = d / np.linalg.norm(d)                    # in-plane, already unit up to rounding
    if loop["shape"] == "circular":
        pieces = [_clip_disk(q, d, np.zeros(2), loop["radius"])]
    else:
        rc = min(loop["corner_radius"], loop["half_m"], loop["half_p"])
        im, ip = loop["half_m"] - rc, loop["half_p"] - rc
        pieces = [_clip_box(q, d, im + rc, ip), _clip_box(q, d, im, ip + rc)]
        if rc > 0:
            pieces += [_clip_disk(q, d, np.array([su * im, sv * ip]), rc) for su in (-1, 1) for sv in (-1, 1)]
    pieces = [pc for pc in pieces if pc is not None]
    if not pieces:
        return None
    return min(pc[0] for pc in pieces), max(pc[1] for pc in pieces)


# =============================================================================
# Writers
# =============================================================================
def _g(x):
    x = 0.0 if abs(x) < 1e-9 else float(x)
    return f"{x:.10g}"


def write_yaml(path, general, loops, header_lines=()):
    """exec_cuda input with the exact key set accepted by input.cpp."""
    lines = [f"# {h}" if h else "#" for h in header_lines]
    if lines:
        lines.append("")
    lines += ["general:", f"  a0: {general['a0']}", f"  Nx: {general['Nx']}", f"  Ny: {general['Ny']}",
              f"  Nz: {general['Nz']}", f"  Loopon_size: {_g(general['Loopon_size'])}"]
    if general.get("Output_format"):
        lines.append(f"  Output_format: {general['Output_format']}")
    lines += ["", "dislocations:"]
    for l in loops:
        lines.append(f"  - Shape: {l['shape']}")
        if l["shape"] == "circular":
            lines.append(f"    Radius: {_g(l['radius'])}")
        elif l["shape"] == "parallelogram":
            lines += [f"    Half_1: {_g(l['half_1'])}", f"    Half_2: {_g(l['half_2'])}",
                      f"    Dir_1: [{', '.join(str(int(x)) for x in l['dir_1'])}]",
                      f"    Dir_2: [{', '.join(str(int(x)) for x in l['dir_2'])}]"]
        else:
            lines += [f"    Half_m: {_g(l['half_m'])}", f"    Half_p: {_g(l['half_p'])}",
                      f"    Corner_radius: {_g(l['corner_radius'])}"]
            if l.get("side") is not None:
                lines.append(f"    Side: [{', '.join(str(int(x)) for x in l['side'])}]")
        lines += [f"    Center: [{', '.join(_g(x) for x in l['center'])}]",
                  f"    m: [{', '.join(str(int(x)) for x in l['m'])}]",
                  f"    n: [{', '.join(str(int(x)) for x in l['n'])}]"]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return round(float(x), 8) + 0.0
    return x


def loop_record(bcc, loop, general):
    """Per-loop ground truth for the sidecar."""
    a0 = general["a0"]
    cells = (general["Nx"], general["Ny"], general["Nz"])
    ty, idx = expected_type_and_sys(bcc, loop["m"], loop["n"])
    rec = dict(shape=loop["shape"])
    if loop["shape"] == "circular":
        rec["radius"] = loop["radius"]
    elif loop["shape"] == "parallelogram":
        rec.update(half_1=loop["half_1"], half_2=loop["half_2"], dir_1=list(loop["dir_1"]),
                   dir_2=list(loop["dir_2"]))
    else:
        rec.update(half_m=loop["half_m"], half_p=loop["half_p"], corner_radius=loop["corner_radius"])
        if loop.get("side") is not None:
            rec["side"] = list(loop["side"])
    rec.update(center_requested=list(loop["center"]),
               center_snapped=list(snap_center(loop["center"], loop["n"], a0, cells)),
               m=list(loop["m"]), n=list(loop["n"]), burgers_A=[a0 / 2 * c for c in loop["m"]],
               plane_family=plane_family(loop["n"]), expected_burgers_type=ty, expected_slip_system=idx,
               loopons=count_loopons(loop, general["Loopon_size"]))
    return rec


def write_sidecar(path, bcc, general, loops, expected, extra=None, script=None):
    """JSON ground truth next to the YAML: generator command line, general block,
    per-loop records and the expected DXA/DANCE content of the sample."""
    doc = dict(generator=script or os.path.basename(sys.argv[0]),
               command=" ".join([os.path.basename(sys.argv[0])] + sys.argv[1:]),
               general=general, loops=[loop_record(bcc, l, general) for l in loops], expected=expected)
    if extra:
        doc.update(extra)
    text = json.dumps(_jsonable(doc), indent=2)
    # short numeric vectors on one line
    text = re.sub(r"\[\s*(-?[\d.e+-]+(?:,\s*-?[\d.e+-]+){1,3})\s*\]",
                  lambda mt: "[" + ", ".join(x.strip() for x in mt.group(1).split(",")) + "]", text)
    with open(path, "w") as f:
        f.write(text + "\n")


def general_block(a0, cells, loopon_size, fmt):
    if isinstance(cells, int):
        cells = (cells, cells, cells)
    return dict(a0=a0, Nx=cells[0], Ny=cells[1], Nz=cells[2], Loopon_size=loopon_size, Output_format=fmt)
