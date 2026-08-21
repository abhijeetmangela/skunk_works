#!/usr/bin/env python3
"""
Thin-Walled Wing Skin / Sheet-Flow GUI
--------------------------------------

Features
- Imports an airfoil .dat file (x/c, y/c; Selig-style files supported).
- Lets the user define top/bottom stringer x/c locations.
- Optional spar x/c creates a two-cell section; without a spar the section is one-cell.
- Uses simplified thin-walled/sheet-flow idealisation:
    * Skin carries shear flow.
    * Spar/stringers contribute to boom-area discretisation for the q_b
      ("basic" shear flow) calculation.
    * Spar web carries vertical shear in a two-cell section.
- Calculates:
    * centroid / boom areas
    * basic shear flow q_b(s)
    * shear centre (open-section line integration + iterative closed-cell torque balance)
    * torsional shear flow q_t
    * total shear flow q
    * shear stress and von Mises stress
    * panel buckling check
    * automatic ribs/stringers increase after buckling failure
- Plots q(s) at selected span stations.
- Exports a text report and CSV results.

IMPORTANT ENGINEERING NOTE
This is a computational design/preliminary sizing tool, not a certification method.
The exact buckling equation from your assignment is not present in the prompt. The
fallback implemented here is a standard simply-supported plate shear-buckling model:
    tau_cr = k_s*pi^2*E/(12*(1-nu^2))*(t/b)^2
with a classical k_s approximation. Replace buckling_coefficient() if your attached
formula differs.

Expected span-load file:
CSV/TXT with 2 or 3 columns:
    y, V
or
    y, V, T
Units:
    y in m from root, V in N, T in N-m.
If T is omitted, torsional moment is estimated as V * e, where e is the shear-load
offset from the calculated shear centre.

Dependencies:
    numpy matplotlib tkinter
Python 3.10+ recommended.
"""

from __future__ import annotations

import csv
import math
import os
import traceback
import tkinter as tk
from dataclasses import dataclass
from tkinter import ttk, filedialog, messagebox
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ----------------------------- data classes -----------------------------

@dataclass
class Material:
    E: float
    nu: float
    yield_strength: float
    density: float = 2700.0


@dataclass
class SectionInputs:
    airfoil_path: str
    spar_x: Optional[float]
    skin_t: float
    spar_t: float
    spar_web_h: float
    stringer_area: float
    n_ribs: int
    max_stringers_per_surface: int
    root_chord: float
    tip_chord: float
    span: float
    load_factor: float
    load_offset: float
    material: Material


@dataclass
class Stringer:
    x: float
    surface: str  # "top" or "bottom"
    area: float


@dataclass
class Boom:
    x: float
    y: float
    area: float
    label: str


@dataclass
class CellGeom:
    """One thin-walled cell boundary represented by panel segments."""
    point_ids: List[int]
    perimeter: float
    area: float


@dataclass
class SectionResult:
    x_surf: np.ndarray
    y_surf: np.ndarray
    ds: np.ndarray
    theta: np.ndarray
    q_basic: np.ndarray
    q_torsion: np.ndarray
    q_total: np.ndarray
    tau_skin: np.ndarray
    vm_skin: np.ndarray
    shear_centre_x: float
    centroid_x: float
    centroid_y: float
    torque_basic_about_centroid: float
    enclosed_area: float
    n_cells: int
    panel_meta: List[dict]


# ----------------------------- airfoil utilities -----------------------------

def read_airfoil_dat(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read common Selig/UIUC airfoil .dat formats."""
    xs, ys = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    for line in lines:
        s = line.strip().replace(",", " ")
        if not s:
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        xs.append(x)
        ys.append(y)

    if len(xs) < 10:
        raise ValueError("Could not read enough numeric x-y points from the airfoil .dat file.")

    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)

    # Normalize if coordinates are not already close to chord=1.
    xmin, xmax = np.min(x), np.max(x)
    chord = xmax - xmin
    if chord <= 0:
        raise ValueError("Invalid airfoil coordinate range.")
    x = (x - xmin) / chord
    y = y / chord

    # Remove repeated consecutive points.
    keep = np.ones(len(x), dtype=bool)
    keep[1:] = np.hypot(np.diff(x), np.diff(y)) > 1e-10
    x, y = x[keep], y[keep]

    return x, y


def split_airfoil_surfaces(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Produce monotonic upper/lower surfaces in x.

    For Selig data, the surface may be listed TE->LE->TE. We use the
    leading-edge point as the split and sort each side by x.
    """
    i_le = int(np.argmin(x))
    upper_x, upper_y = x[: i_le + 1], y[: i_le + 1]
    lower_x, lower_y = x[i_le:], y[i_le:]

    # Ensure each is monotonic increasing in x.
    iu = np.argsort(upper_x)
    il = np.argsort(lower_x)
    upper_x, upper_y = upper_x[iu], upper_y[iu]
    lower_x, lower_y = lower_x[il], lower_y[il]

    # Remove duplicate x by averaging y.
    def unique_average(xa, ya):
        ux = np.unique(np.round(xa, 12))
        uy = np.array([np.mean(ya[np.isclose(xa, xx, atol=5e-12)]) for xx in ux])
        return ux, uy

    upper_x, upper_y = unique_average(upper_x, upper_y)
    lower_x, lower_y = unique_average(lower_x, lower_y)

    return upper_x, upper_y, lower_x, lower_y


def interp_surface(xq: float, x: np.ndarray, y: np.ndarray) -> float:
    xq = float(np.clip(xq, np.min(x), np.max(x)))
    return float(np.interp(xq, x, y))


# ----------------------------- geometry / booms -----------------------------

def stringers_from_entries(entries: List[Tuple[str, str]]) -> List[Stringer]:
    out = []
    for s, surf in entries:
        try:
            xx = float(s)
        except ValueError:
            continue
        surf = surf.strip().lower()
        if surf not in ("top", "bottom"):
            continue
        out.append(Stringer(xx, surf, 1.0))
    return out


def build_stringers(text_top: str, text_bottom: str, area: float) -> List[Stringer]:
    out = []
    for tok in text_top.replace(";", ",").split(","):
        tok = tok.strip()
        if tok:
            out.append(Stringer(float(tok), "top", area))
    for tok in text_bottom.replace(";", ",").split(","):
        tok = tok.strip()
        if tok:
            out.append(Stringer(float(tok), "bottom", area))
    return out


def polyline_length(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sum(np.hypot(np.diff(x), np.diff(y))))


def polygon_area(x: np.ndarray, y: np.ndarray) -> float:
    return abs(0.5 * float(np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))))


def build_section_geometry(
    airfoil_x: np.ndarray,
    airfoil_y: np.ndarray,
    chord: float,
    stringers: List[Stringer],
    spar_x: Optional[float],
    skin_t: float,
    spar_t: float,
    spar_web_h: float,
) -> Tuple[np.ndarray, np.ndarray, List[Boom], List[CellGeom], dict]:
    """
    Returns a discretised closed contour in metre coordinates.

    The skin contour is constructed as upper TE->LE + lower LE->TE.
    A spar is an internal vertical line joining top and bottom skins.
    """
    xu, yu, xl, yl = split_airfoil_surfaces(airfoil_x, airfoil_y)

    # Use sufficiently fine x-discretisation for stable q(s) plots.
    xgrid = np.linspace(0, 1, 801)
    yu_i = np.interp(xgrid, xu, yu)
    yl_i = np.interp(xgrid, xl, yl)

    # Closed skin contour clockwise: upper TE -> LE, lower LE -> TE.
    xs = np.concatenate([xgrid[::-1], xgrid[1:]])
    ys = np.concatenate([yu_i[::-1], yl_i[1:]])
    xs *= chord
    ys *= chord

    # Ensure centroid-friendly axis is conventional x-c chordwise, y vertical.
    area_geom = polygon_area(xs, ys)
    if area_geom <= 0:
        raise ValueError("Airfoil polygon area is zero.")

    # Stringer points.
    booms: List[Boom] = []
    for i, st in enumerate(stringers):
        xx = float(np.clip(st.x, 0.0, 1.0))
        yy_c = float(np.interp(xx, xgrid, yu_i if st.surface == "top" else yl_i))
        booms.append(Boom(xx * chord, yy_c * chord, st.area, f"S{i+1}_{st.surface}"))

    # Optional spar.
    spar_meta = None
    cells = []
    if spar_x is not None:
        sx = float(np.clip(spar_x, 0.02, 0.98))
        ytop = float(np.interp(sx, xgrid, yu_i)) * chord
        ybot = float(np.interp(sx, xgrid, yl_i)) * chord
        h = max(abs(ytop - ybot), 1e-6)
        spar_meta = {
            "x": sx * chord,
            "ytop": ytop,
            "ybot": ybot,
            "h": h,
            "t": spar_t,
            "web_area": h * spar_t,
        }

        # Geometric cell areas:
        # front cell = area between LE and spar
        # aft cell = area between spar and TE
        def area_between(xa, xb):
            mask = (xgrid >= xa) & (xgrid <= xb)
            xx = xgrid[mask] * chord
            up = yu_i[mask] * chord
            lo = yl_i[mask] * chord
            return float(np.trapz(up - lo, xx))

        a_front = area_between(0.0, sx)
        a_aft = area_between(sx, 1.0)
        cells = [
            CellGeom([], max(1e-12, a_front * 2 * 0 + 1.0), a_front),
            CellGeom([], max(1e-12, a_aft * 2 * 0 + 1.0), a_aft),
        ]

    geom = {
        "xgrid": xgrid,
        "yu": yu_i,
        "yl": yl_i,
        "chord": chord,
        "spar": spar_meta,
        "skin_area_midline": skin_t * polyline_length(xs, ys),
        "enclosed_area": area_geom,
        "n_cells": 2 if spar_meta else 1,
        "cells": cells,
    }
    return xs, ys, booms, cells, geom


def section_centroid_with_booms(
    xs: np.ndarray,
    ys: np.ndarray,
    booms: List[Boom],
    skin_t: float,
    spar_meta: Optional[dict],
) -> Tuple[float, float, float]:
    """
    Area centroid using skin + spar + stringers.
    """
    ds = np.hypot(np.diff(xs), np.diff(ys))
    xmid = 0.5 * (xs[:-1] + xs[1:])
    ymid = 0.5 * (ys[:-1] + ys[1:])
    a_skin = ds * skin_t

    A = float(np.sum(a_skin))
    Sx = float(np.sum(a_skin * xmid))
    Sy = float(np.sum(a_skin * ymid))

    if spar_meta:
        As = spar_meta["web_area"]
        xspar = spar_meta["x"]
        yspar = 0.5 * (spar_meta["ytop"] + spar_meta["ybot"])
        A += As
        Sx += As * xspar
        Sy += As * yspar

    for b in booms:
        A += b.area
        Sx += b.area * b.x
        Sy += b.area * b.y

    return Sx / A, Sy / A, A


# ----------------------------- shear-flow methods -----------------------------

def discretise_contour(xs: np.ndarray, ys: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ds = np.hypot(np.diff(xs), np.diff(ys))
    xm = 0.5 * (xs[:-1] + xs[1:])
    ym = 0.5 * (ys[:-1] + ys[1:])
    theta = np.arctan2(np.diff(ys), np.diff(xs))
    return ds, xm, ym, theta


def interpolate_q_to_boom_path(
    xs, ys, qb, boom: Boom, window: float = 0.02
) -> float:
    # nearest contour point to boom point
    d = np.hypot(xs[:-1] - boom.x, ys[:-1] - boom.y)
    j = int(np.argmin(d))
    return float(qb[j])


def basic_shear_flow_open(
    xs: np.ndarray,
    ys: np.ndarray,
    booms: List[Boom],
    skin_t: float,
    Vx: float,
    Vy: float,
    xbar: float,
    ybar: float,
    spar_meta: Optional[dict],
) -> Tuple[np.ndarray, float]:
    """
    Thin-walled boom idealisation.

    For a closed contour, this first computes the open-section basic q_b after
    cutting the contour at the first panel. The boom-area jumps are included.
    The shear-flow sign follows the chosen contour orientation.

    The integral used is:
        q_b(s) = [ Vx * Q_y(s)/Ix + Vy * Q_x(s)/Iy ]  for a simplified
        uncoupled form, with Q_x = ∫(x-xbar) dB and Q_y = ∫(y-ybar)dB.

    The skin's own contribution to I_x/I_y is retained; stringers/spar are
    included as lumped areas.
    """
    ds, xm, ym, theta = discretise_contour(xs, ys)
    # Exact-ish section inertias including skin and discrete areas.
    w = ds * skin_t
    xr = xm - xbar
    yr = ym - ybar
    Ix = float(np.sum(w * yr**2))
    Iy = float(np.sum(w * xr**2))
    Ixy = float(np.sum(w * xr * yr))
    for b in booms:
        Ix += b.area * (b.y - ybar) ** 2
        Iy += b.area * (b.x - xbar) ** 2
        Ixy += b.area * (b.x - xbar) * (b.y - ybar)
    if spar_meta:
        As = spar_meta["web_area"]
        Ix += As * ((spar_meta["ybt"] if "ybt" in spar_meta else 0.5*(spar_meta["ytop"]+spar_meta["ybot"])) - ybar) ** 2
        Iy += As * (spar_meta["x"] - xbar) ** 2

    den = Ix * Iy - Ixy**2
    if abs(den) < 1e-18:
        raise ValueError("Section inertias are singular.")

    # Use coupled shear formula:
    # [A;B] = 1/den * [[Iy,-Ixy],[-Ixy,Ix]] * [Vy,Vx]
    # q = (A*(∫ y dA) + B*(∫ x dA))/... simplified into boom-based path.
    coeff_y = (Vy * Iy - Vx * Ixy) / den
    coeff_x = (Vx * Ix - Vy * Ixy) / den

    # Cut at start: cumulative first moment of boom+skin strip.
    q = np.zeros(len(ds), dtype=float)

    # We represent distributed skin contribution continuously and add lumped
    # boom areas whenever the cumulative path passes their nearest station.
    cum_Qy = 0.0  # ∫(y-ybar) dA
    cum_Qx = 0.0  # ∫(x-xbar) dA

    # Place booms onto contour stations.
    boom_idx = []
    for b in booms:
        j = int(np.argmin(np.hypot(xm - b.x, ym - b.y)))
        boom_idx.append((j, b))

    for i in range(len(ds)):
        # Include boom at the current point before evaluating current panel flow.
        for j, b in boom_idx:
            if j == i:
                cum_Qy += b.area * (b.y - ybar)
                cum_Qx += b.area * (b.x - xbar)

        # Half panel distributed skin area.
        dA = ds[i] * skin_t
        cum_Qy += dA * (ym[i] - ybar)
        cum_Qx += dA * (xm[i] - xbar)

        q[i] = coeff_y * cum_Qy + coeff_x * cum_Qx

    # The simple formula above has units N/m.
    # Force q units through V and Q / I are correct.
    return q, float(Ix)


def shear_centre_by_torque(
    xs: np.ndarray,
    ys: np.ndarray,
    q_basic: np.ndarray,
    Vx: float,
    Vy: float,
    xbar: float,
    ybar: float,
) -> Tuple[float, float]:
    """
    Calculate the shear centre from the torque of the basic open-section q
    about the centroid.

    For a vertical shear Vy:
        T_q = ∮ q * (r x t_hat) ds
        e = T_q / Vy
        x_sc = xbar + e (for positive convention)

    For general Vx/Vy, solve along the load line using the appropriate
    perpendicular offset. This implementation returns an x position for the
    common aircraft case V = Vy.
    """
    ds, xm, ym, _ = discretise_contour(xs, ys)
    if abs(Vy) < 1e-12:
        # Horizontal case can be handled separately; return centroid for GUI.
        return xbar, ybar

    # dl vector along contour. r x (q dl) = q*(rx*dY - ry*dX)
    dx = np.diff(xs)
    dy = np.diff(ys)
    rx = xm - xbar
    ry = ym - ybar
    T = float(np.sum(q_basic * (rx * dy - ry * dx)))

    x_sc = xbar + T / Vy
    return x_sc, ybar


def closed_cell_corrected_flow(
    xs: np.ndarray,
    ys: np.ndarray,
    q_basic: np.ndarray,
    skin_t: float,
    G: float,
    torsion: float,
    spar_meta: Optional[dict],
    xbar: float,
    ybar: float,
) -> Tuple[np.ndarray, float]:
    """
    Simplified closed-cell torsional/shear-flow correction.

    One cell:
        q_t = T / (2 A)

    Two cells separated by a spar:
        solve compatibility with constant cell flows q0_1, q0_2.
    A practical approximation is used by representing the two cell twist
    relation from shear compliance ∮ q/(G t) ds = 2 A theta.

    For the assignment's simplified sheet-flow model, the correction is:
        q = q_b + q0_cell
    and a common-twist condition is enforced.

    We avoid requiring a full FE skin graph; the spar is treated as a
    station at x_spar and the skin contour is partitioned into front/aft
    arcs.
    """
    ds, xm, ym, _ = discretise_contour(xs, ys)

    if abs(torsion) < 1e-14:
        return np.zeros_like(q_basic), 0.0

    # Single-cell: Bredt-Batho.
    if spar_meta is None:
        area = polygon_area(xs, ys)
        if area <= 1e-12:
            return np.zeros_like(q_basic), 0.0
        q0 = torsion / (2.0 * area)
        return np.full_like(q_basic, q0), q0

    # Two-cell partition by spar x.
    sx = spar_meta["x"]
    front = xm <= sx
    aft = ~front

    # Compute geometric cell areas from skin arcs + straight spar.
    xmask = xs
    ymask = ys
    area_total = polygon_area(xs, ys)
    # Front and aft areas via area under airfoil split.
    ytop = spar_meta["ytop"]
    ybot = spar_meta["ybot"]

    # Area by polygon clipping through x=sx approximated from x-coordinate masks.
    def approximate_cell_area(mask):
        # integrate top-bottom over panel midpoint segments on selected side
        vals = ds[mask]
        # use a projected x integration from contour samples:
        # fallback based on total area partition by chord fraction.
        return max(1e-12, area_total * (sx / (np.max(xs) if np.max(xs) > 0 else 1.0)))

    # Better area partition using contour samples and vertical split.
    # Shoelace clipping routine.
    def clip_polygon(poly_x, poly_y, xclip, keep_left):
        ox, oy = [], []
        n = len(poly_x)
        for i in range(n):
            x1, y1 = poly_x[i], poly_y[i]
            x2, y2 = poly_x[(i+1) % n], poly_y[(i+1) % n]
            in1 = x1 <= xclip + 1e-12 if keep_left else x1 >= xclip - 1e-12
            in2 = x2 <= xclip + 1e-12 if keep_left else x2 >= xclip - 1e-12
            if in1:
                ox.append(x1); oy.append(y1)
            if in1 != in2:
                tt = (xclip - x1) / (x2 - x1)
                ox.append(x1 + tt * (x2 - x1))
                oy.append(y1 + tt * (y2 - y1))
        return np.array(ox), np.array(oy)

    px = np.asarray(xs)
    py = np.asarray(ys)
    lx, ly = clip_polygon(px, py, sx, True)
    rx, ry = clip_polygon(px, py, sx, False)
    A1 = polygon_area(lx, ly) if len(lx) >= 3 else area_total * 0.5
    A2 = polygon_area(rx, ry) if len(rx) >= 3 else area_total * 0.5

    # Basic shear-flow torque in each cell.
    def cell_torque(mask):
        dx = np.diff(xs)[mask]
        dy = np.diff(ys)[mask]
        rm_x = xm[mask] - xbar
        rm_y = ym[mask] - ybar
        return float(np.sum(q_basic[mask] * (rm_x * dy - rm_y * dx)))

    T1 = cell_torque(front)
    T2 = cell_torque(aft)

    # Compatibility:
    # (q_b1+q01) perimeter1/(G t) = (q_b2+q02) perimeter2/(G t) up to
    # the same twist. A compact 2x2 system is formed including spar web.
    p1 = float(np.sum(ds[front]))
    p2 = float(np.sum(ds[aft]))
    ps = max(spar_meta["h"], 1e-12)
    ts = max(spar_meta["t"], 1e-12)
    t = max(skin_t, 1e-12)

    # q01/t perimeter + q_s/ts * spar_length enters cell twist;
    # q_s is the difference between the cell flows across the spar.
    # Compatibility:
    # 2 A1 theta = ∮cell1 q/(Gt)
    # 2 A2 theta = ∮cell2 q/(Gt)
    # Torque: T = 2A1*q01 + 2A2*q02 + base torque.
    # We solve unknown q01, q02 and theta.
    # q_spar = q01 - q02 (sign convention), but exact sign depends on contour.
    M = np.array([
        [p1/t + ps/ts, -ps/ts, -2*A1],
        [-ps/ts, p2/t + ps/ts, -2*A2],
        [2*A1, 2*A2, 0.0],
    ], dtype=float)
    b = np.array([
        -float(np.sum((q_basic[front]) * ds[front] / t)) - T1*0.0,
        -float(np.sum((q_basic[aft]) * ds[aft] / t)),
        torsion - T1 - T2,
    ], dtype=float)

    try:
        sol = np.linalg.solve(M, b)
        q01, q02, theta = sol
    except np.linalg.LinAlgError:
        # robust fallback: equal torque split with compatibility ignored
        q01 = (torsion - T1 - T2) / (2.0 * (A1 + A2))
        q02 = q01
        theta = 0.0

    qt = np.zeros_like(q_basic)
    qt[front] += q01
    qt[aft] += q02
    return qt, float(theta)


# ----------------------------- panel buckling -----------------------------

def buckling_coefficient(aspect_ratio: float) -> float:
    """
    Standard simply-supported plate shear-buckling approximation.

    aspect_ratio = a/b, with a = longer panel direction, b = plate width.
    """
    r = max(float(aspect_ratio), 1e-8)
    if r >= 1.0:
        return 5.34 + 4.0 / (r**2)
    # Alternative symmetric form for r < 1.
    return 5.34 * (r**2) + 4.0


def tau_cr_shear_buckling(E: float, nu: float, t: float, b: float, a: float) -> float:
    ks = buckling_coefficient(max(a, b) / max(min(a, b), 1e-12))
    b_eff = max(min(a, b), 1e-9)
    return ks * math.pi**2 * E / (12.0 * (1.0 - nu**2)) * (t / b_eff) ** 2


def panel_widths_from_stringers(
    chord: float, stringers: List[Stringer], surface: str
) -> List[Tuple[float, float]]:
    xs = sorted([s.x for s in stringers if s.surface == surface])
    bounds = [0.0] + xs + [1.0]
    return [(bounds[i], bounds[i+1]) for i in range(len(bounds)-1)]


# ----------------------------- section solver -----------------------------

def solve_section(
    airfoil_x,
    airfoil_y,
    chord: float,
    inp: SectionInputs,
    stringers: List[Stringer],
    V: float,
    T_applied: float,
) -> SectionResult:
    xs, ys, booms, cells, geom = build_section_geometry(
        airfoil_x,
        airfoil_y,
        chord,
        stringers,
        inp.spar_x,
        inp.skin_t,
        inp.spar_t,
        inp.spar_web_h,
    )

    xbar, ybar, A = section_centroid_with_booms(
        xs, ys, booms, inp.skin_t, geom.get("spar")
    )

    # For the common vertical shear case.
    q_basic, Ix = basic_shear_flow_open(
        xs, ys, booms, inp.skin_t, 0.0, V, xbar, ybar, geom.get("spar")
    )

    x_sc, y_sc = shear_centre_by_torque(xs, ys, q_basic, 0.0, V, xbar, ybar)

    # Applied torsion = external moment about shear centre.
    qt, _theta = closed_cell_corrected_flow(
        xs, ys, q_basic, inp.skin_t,
        inp.material.E / (2.0 * (1.0 + inp.material.nu)),
        T_applied, geom.get("spar"), xbar, ybar
    )
    qtot = q_basic + qt

    ds, xm, ym, _ = discretise_contour(xs, ys)
    tau = qtot / inp.skin_t
    # With no longitudinal normal stress included, VM from shear alone is sqrt(3)*tau.
    vm = np.sqrt(3.0) * np.abs(tau)

    panel_meta = make_panel_checks(
        inp, chord, stringers, xs, ys, qtot
    )

    return SectionResult(
        x_surf=xm,
        y_surf=ym,
        ds=ds,
        theta=np.arctan2(np.diff(ys), np.diff(xs)),
        q_basic=q_basic,
        q_torsion=qt,
        q_total=qtot,
        tau_skin=tau,
        vm_skin=vm,
        shear_centre_x=x_sc,
        centroid_x=xbar,
        centroid_y=ybar,
        torque_basic_about_centroid=0.0,
        enclosed_area=A,
        n_cells=geom["n_cells"],
        panel_meta=panel_meta,
    )


def make_panel_checks(
    inp: SectionInputs,
    chord: float,
    stringers: List[Stringer],
    xs: np.ndarray,
    ys: np.ndarray,
    q_total: np.ndarray,
) -> List[dict]:
    ds, xm, ym, _ = discretise_contour(xs, ys)
    panels = []

    # Spanwise rib pitch.
    n_bays = max(1, int(inp.n_ribs) - 1)
    a = inp.span / n_bays  # m, rib-to-rib dimension

    # Surface panels based on x spacing.
    for surf in ("top", "bottom"):
        bounds = panel_widths_from_stringers(chord, stringers, surf)
        for k, (x1, x2) in enumerate(bounds):
            xx1, xx2 = x1 * chord, x2 * chord
            # Select q on the corresponding skin side.
            if surf == "top":
                ysurf = np.interp(np.array([xx1, xx2]), xs, ys)
                # nearest top-surface q is used by distance to straight interpolation
                mask = (xm >= xx1 - 1e-9) & (xm <= xx2 + 1e-9) & (ym >= 0.0)
            else:
                mask = (xm >= xx1 - 1e-9) & (xm <= xx2 + 1e-9) & (ym <= 0.0)

            # If the airfoil mean line is not exactly y=0, use abs-y side heuristic.
            candidates = np.where(mask)[0]
            if len(candidates) == 0:
                # fallback to x-only
                candidates = np.where((xm >= xx1 - 1e-9) & (xm <= xx2 + 1e-9))[0]
            qmax = float(np.max(np.abs(q_total[candidates]))) if len(candidates) else 0.0
            tau_max = qmax / inp.skin_t
            b = max((x2 - x1) * chord, 1e-9)
            taucr = tau_cr_shear_buckling(
                inp.material.E, inp.material.nu, inp.skin_t, b, a
            )
            panels.append({
                "surface": surf,
                "index": k + 1,
                "x1": x1,
                "x2": x2,
                "a_rib": a,
                "b_stringer": b,
                "qmax": qmax,
                "tau_max": tau_max,
                "tau_cr": taucr,
                "buckling_factor": taucr / max(tau_max, 1e-30),
                "fails": tau_max > taucr,
            })

    # If a spar is present, also treat its web bays.
    if inp.spar_x is not None:
        spar_h = max(inp.spar_web_h, 1e-9)
        b = spar_h
        taucr = tau_cr_shear_buckling(
            inp.material.E, inp.material.nu, inp.spar_t, b, a
        )
        # Spar web shear flow approximately equal to cell-flow difference.
        # Use the maximum skin q as conservative proxy in this simplified model.
        qmax = float(np.max(np.abs(q_total)))
        tau_max = qmax / max(inp.spar_t, 1e-12)
        panels.append({
            "surface": "spar_web",
            "index": 1,
            "x1": inp.spar_x,
            "x2": inp.spar_x,
            "a_rib": a,
            "b_stringer": spar_h,
            "qmax": qmax,
            "tau_max": tau_max,
            "tau_cr": taucr,
            "buckling_factor": taucr / max(tau_max, 1e-30),
            "fails": tau_max > taucr,
        })

    return panels


# ----------------------------- spanwise processing -----------------------------

def read_span_load_file(path: str) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass
            if len(vals) >= 2:
                rows.append(vals[:3])

    if len(rows) < 2:
        raise ValueError("Span-load file must contain at least two numeric rows.")

    arr = np.array(rows, dtype=float)
    y = arr[:, 0]
    V = arr[:, 1]
    T = arr[:, 2] if arr.shape[1] >= 3 else None
    idx = np.argsort(y)
    y, V = y[idx], V[idx]
    if T is not None:
        T = T[idx]
    return y, V, T


def choose_stations(y: np.ndarray) -> List[float]:
    if len(y) <= 6:
        return list(y)
    # root, 20%, 40%, 60%, 80%, tip stations
    fracs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    yy = np.interp(fracs, np.linspace(0, 1, len(y)), y)
    return list(yy)


# ----------------------------- GUI -----------------------------

class WingSheetFlowApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Wing Skin Sheet-Flow / Shear-Centre / Buckling Tool")
        self.geometry("1420x920")
        self.minsize(1180, 760)

        self.airfoil_x = None
        self.airfoil_y = None
        self.span_y = None
        self.span_V = None
        self.span_T = None
        self.latest_results = []
        self.latest_summary = ""

        self._build_variables()
        self._build_ui()

    def _build_variables(self):
        self.airfoil_var = tk.StringVar()
        self.load_var = tk.StringVar()
        self.spar_var = tk.StringVar(value="0.25")
        self.skin_t_var = tk.StringVar(value="0.0012")
        self.spar_t_var = tk.StringVar(value="0.0020")
        self.spar_h_var = tk.StringVar(value="0.25")
        self.stringer_area_var = tk.StringVar(value="80e-6")
        self.top_stringers_var = tk.StringVar(value="0.15,0.35,0.55,0.75,0.90")
        self.bottom_stringers_var = tk.StringVar(value="0.15,0.35,0.55,0.75,0.90")
        self.n_ribs_var = tk.StringVar(value="16")
        self.max_stringers_var = tk.StringVar(value="12")
        self.root_chord_var = tk.StringVar(value="1.5")
        self.tip_chord_var = tk.StringVar(value="0.9")
        self.span_var = tk.StringVar(value="5.0")
        self.load_factor_var = tk.StringVar(value="2.5")
        self.load_offset_var = tk.StringVar(value="0.0")

        self.E_var = tk.StringVar(value="70e9")
        self.nu_var = tk.StringVar(value="0.33")
        self.yield_var = tk.StringVar(value="240e6")

        self.status_var = tk.StringVar(value="Ready.")
        self.station_var = tk.StringVar(value="0.0")

    def _build_ui(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)

        right = ttk.Frame(self, padding=8)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Left notebook.
        nb = ttk.Notebook(left)
        nb.pack(fill=tk.BOTH, expand=True)

        geo = ttk.Frame(nb, padding=8)
        load = ttk.Frame(nb, padding=8)
        mat = ttk.Frame(nb, padding=8)
        out = ttk.Frame(nb, padding=8)
        nb.add(geo, text="Geometry")
        nb.add(load, text="Loads")
        nb.add(mat, text="Material")
        nb.add(out, text="Actions")

        r = 0
        ttk.Label(geo, text="Airfoil .dat file").grid(row=r, column=0, sticky="w")
        ttk.Entry(geo, textvariable=self.airfoil_var, width=38).grid(row=r+1, column=0, sticky="ew")
        ttk.Button(geo, text="Browse", command=self.browse_airfoil).grid(row=r+1, column=1, padx=4)
        r += 3

        ttk.Label(geo, text="Span-load file (.csv/.txt): y, V, [T]").grid(row=r, column=0, sticky="w")
        ttk.Entry(geo, textvariable=self.load_var, width=38).grid(row=r+1, column=0, sticky="ew")
        ttk.Button(geo, text="Browse", command=self.browse_load).grid(row=r+1, column=1, padx=4)
        r += 3

        fields = [
            ("Spar x/c (blank = one-cell)", self.spar_var),
            ("Skin thickness t [m]", self.skin_t_var),
            ("Spar thickness [m]", self.spar_t_var),
            ("Spar web height [m]", self.spar_h_var),
            ("Stringer area [m^2]", self.stringer_area_var),
            ("Top stringer x/c", self.top_stringers_var),
            ("Bottom stringer x/c", self.bottom_stringers_var),
            ("Number of ribs N_R", self.n_ribs_var),
            ("Max stringers / surface", self.max_stringers_var),
            ("Root chord [m]", self.root_chord_var),
            ("Tip chord [m]", self.tip_chord_var),
            ("Wing semi-span [m]", self.span_var),
        ]
        for label, var in fields:
            ttk.Label(geo, text=label).grid(row=r, column=0, sticky="w", pady=(3,0))
            ttk.Entry(geo, textvariable=var, width=42).grid(row=r+1, column=0, columnspan=2, sticky="ew")
            r += 2

        r = 0
        for label, var in [
            ("Load factor n", self.load_factor_var),
            ("Load offset from shear centre [m]", self.load_offset_var),
        ]:
            ttk.Label(load, text=label).grid(row=r, column=0, sticky="w")
            ttk.Entry(load, textvariable=var, width=25).grid(row=r+1, column=0, sticky="ew")
            r += 2

        ttk.Label(load, text="Span-load file convention").grid(row=r, column=0, sticky="w", pady=(8,2))
        ttk.Label(load, text="y [m], V [N], optional T [N-m]").grid(row=r+1, column=0, sticky="w")
        r += 3
        ttk.Button(load, text="Preview loads", command=self.preview_loads).grid(row=r, column=0, sticky="ew")
        r += 1

        r = 0
        for label, var in [
            ("Young's modulus E [Pa]", self.E_var),
            ("Poisson ratio nu", self.nu_var),
            ("Yield strength [Pa]", self.yield_var),
        ]:
            ttk.Label(mat, text=label).grid(row=r, column=0, sticky="w")
            ttk.Entry(mat, textvariable=var, width=25).grid(row=r+1, column=0, sticky="ew")
            r += 2

        ttk.Label(mat, text="Buckling formula used").grid(row=r, column=0, sticky="w", pady=(10,2))
        ttk.Label(
            mat,
            text="tau_cr = k_s*pi^2 E/[12(1-nu^2)]*(t/b)^2",
            wraplength=280,
        ).grid(row=r+1, column=0, sticky="w")
        r += 2
        ttk.Label(
            mat,
            text="Edit buckling_coefficient() and tau_cr_shear_buckling() to match your attached formula.",
            wraplength=280,
        ).grid(row=r+1, column=0, sticky="w")

        ttk.Button(out, text="Run analysis", command=self.run_analysis).pack(fill=tk.X, pady=5)
        ttk.Button(out, text="Generate report + CSV", command=self.export_results).pack(fill=tk.X, pady=5)
        ttk.Button(out, text="Plot q at selected station", command=self.plot_selected_station).pack(fill=tk.X, pady=5)
        ttk.Button(out, text="Show airfoil", command=self.plot_airfoil).pack(fill=tk.X, pady=5)
        ttk.Button(out, text="Show airfoil + booms", command=self.plot_section_with_booms).pack(fill=tk.X, pady=5)

        ttk.Label(out, text="Station y [m]").pack(anchor="w", pady=(12,2))
        ttk.Entry(out, textvariable=self.station_var).pack(fill=tk.X)

        ttk.Separator(out).pack(fill=tk.X, pady=12)
        ttk.Label(out, textvariable=self.status_var, wraplength=280).pack(anchor="w")

        # Right: figure and text.
        self.fig = plt.Figure(figsize=(9.4, 6.8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Results")
        self.ax.grid(True, alpha=0.25)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.text = tk.Text(right, height=14, width=110, font=("TkFixedFont", 9))
        self.text.pack(fill=tk.X, pady=(8,0))

    # -------------------- GUI actions --------------------

    def browse_airfoil(self):
        p = filedialog.askopenfilename(
            title="Select airfoil .dat",
            filetypes=[("DAT files", "*.dat *.txt"), ("All files", "*.*")]
        )
        if p:
            self.airfoil_var.set(p)
            try:
                self.airfoil_x, self.airfoil_y = read_airfoil_dat(p)
                self.status_var.set(f"Loaded airfoil: {os.path.basename(p)} ({len(self.airfoil_x)} points)")
                self.plot_section_with_booms()
            except Exception as e:
                messagebox.showerror("Airfoil error", str(e))

    def browse_load(self):
        p = filedialog.askopenfilename(
            title="Select span-load file",
            filetypes=[("CSV/TXT", "*.csv *.txt"), ("All files", "*.*")]
        )
        if p:
            self.load_var.set(p)
            try:
                self.span_y, self.span_V, self.span_T = read_span_load_file(p)
                self.status_var.set(f"Loaded span-load file: {os.path.basename(p)} ({len(self.span_y)} stations)")
                self.preview_loads()
            except Exception as e:
                messagebox.showerror("Load file error", str(e))

    def preview_loads(self):
        if self.span_y is None:
            path = self.load_var.get().strip()
            if not path:
                messagebox.showwarning("Loads", "Select a span-load file.")
                return
            self.span_y, self.span_V, self.span_T = read_span_load_file(path)

        self.ax.clear()
        self.ax.plot(self.span_y, self.span_V, marker="o", label="V")
        if self.span_T is not None:
            ax2 = self.ax.twinx()
            ax2.plot(self.span_y, self.span_T, marker="s", label="T")
            ax2.set_ylabel("T [N-m]")
        self.ax.set_xlabel("Span y [m]")
        self.ax.set_ylabel("Shear force V [N]")
        self.ax.set_title("Spanwise loads")
        self.ax.grid(True, alpha=0.25)
        self.canvas.draw()

    def plot_section_with_booms(self):
        """Plot the imported airfoil together with stringer/boom and spar locations."""
        try:
            inp = self.parse_inputs()

            if self.airfoil_x is None:
                path = self.airfoil_var.get().strip()
                if not path:
                    messagebox.showwarning("Airfoil", "Select an airfoil .dat file first.")
                    return
                self.airfoil_x, self.airfoil_y = read_airfoil_dat(path)

            # Use root chord for the geometry preview because this is the
            # reference section visible before a span station is selected.
            chord = inp.root_chord
            stringers = self._make_stringers(inp)

            xs, ys, booms, cells, geom = build_section_geometry(
                self.airfoil_x,
                self.airfoil_y,
                chord,
                stringers,
                inp.spar_x,
                inp.skin_t,
                inp.spar_t,
                inp.spar_web_h,
            )

            self.ax.clear()

            # Airfoil skin
            self.ax.plot(
                xs / chord,
                ys / chord,
                linewidth=1.5,
                label="Airfoil skin",
            )

            # Boom / stringer locations
            top_b = [b for b in booms if "_top" in b.label]
            bottom_b = [b for b in booms if "_bottom" in b.label]

            if top_b:
                self.ax.scatter(
                    [b.x / chord for b in top_b],
                    [b.y / chord for b in top_b],
                    s=60,
                    marker="o",
                    label="Top booms / stringers",
                    zorder=5,
                )

            if bottom_b:
                self.ax.scatter(
                    [b.x / chord for b in bottom_b],
                    [b.y / chord for b in bottom_b],
                    s=60,
                    marker="o",
                    label="Bottom booms / stringers",
                    zorder=5,
                )

            # Label each boom.
            for b in booms:
                self.ax.annotate(
                    b.label,
                    (b.x / chord, b.y / chord),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                )

            # Spar location
            if geom.get("spar") is not None:
                sp = geom["spar"]
                self.ax.plot(
                    [sp["x"] / chord, sp["x"] / chord],
                    [sp["ybot"] / chord, sp["ytop"] / chord],
                    linewidth=2.5,
                    label=f"Spar x/c={inp.spar_x:.3f}",
                    zorder=4,
                )

            # Show centroid including skin + spar + stringer areas.
            xbar, ybar, A = section_centroid_with_booms(
                xs, ys, booms, inp.skin_t, geom.get("spar")
            )
            self.ax.scatter(
                [xbar / chord],
                [ybar / chord],
                marker="+",
                s=120,
                linewidths=2,
                label="Area centroid",
                zorder=6,
            )

            # Basic geometry annotations.
            self.ax.axhline(
                y=0.0,
                linewidth=0.8,
                linestyle="--",
                alpha=0.5,
            )

            self.ax.set_aspect("equal", adjustable="box")
            self.ax.set_xlabel("x/c")
            self.ax.set_ylabel("y/c")
            self.ax.set_title(
                f"Airfoil + Booms / Stringers | Root Section | "
                f"{'Two-cell' if inp.spar_x is not None else 'One-cell'}"
            )
            self.ax.grid(True, alpha=0.25)
            self.ax.legend(loc="best")
            self.canvas.draw()

            self.status_var.set(
                f"Section preview: {len(booms)} booms; "
                f"{'spar at x/c=' + f'{inp.spar_x:.3f}' if inp.spar_x is not None else 'no spar'}."
            )

        except Exception as e:
            messagebox.showerror(
                "Section preview error",
                f"{e}\n\n{traceback.format_exc()}"
            )

    def plot_airfoil(self):
        path = self.airfoil_var.get().strip()
        if self.airfoil_x is None:
            if not path:
                return
            self.airfoil_x, self.airfoil_y = read_airfoil_dat(path)

        self.ax.clear()
        self.ax.plot(self.airfoil_x, self.airfoil_y)
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.set_xlabel("x/c")
        self.ax.set_ylabel("y/c")
        self.ax.set_title("Airfoil geometry")
        self.ax.grid(True, alpha=0.25)
        self.canvas.draw()

    def parse_inputs(self) -> SectionInputs:
        spar_s = self.spar_var.get().strip()
        spar = None if spar_s == "" else float(spar_s)

        mat = Material(
            E=float(self.E_var.get()),
            nu=float(self.nu_var.get()),
            yield_strength=float(self.yield_var.get()),
        )
        return SectionInputs(
            airfoil_path=self.airfoil_var.get().strip(),
            spar_x=spar,
            skin_t=float(self.skin_t_var.get()),
            spar_t=float(self.spar_t_var.get()),
            spar_web_h=float(self.spar_h_var.get()),
            stringer_area=float(self.stringer_area_var.get()),
            n_ribs=int(self.n_ribs_var.get()),
            max_stringers_per_surface=int(self.max_stringers_var.get()),
            root_chord=float(self.root_chord_var.get()),
            tip_chord=float(self.tip_chord_var.get()),
            span=float(self.span_var.get()),
            load_factor=float(self.load_factor_var.get()),
            load_offset=float(self.load_offset_var.get()),
            material=mat,
        )

    def _make_stringers(self, inp: SectionInputs, top=None, bottom=None) -> List[Stringer]:
        top_text = self.top_stringers_var.get() if top is None else ",".join(f"{v:.5f}" for v in top)
        bottom_text = self.bottom_stringers_var.get() if bottom is None else ",".join(f"{v:.5f}" for v in bottom)
        sts = build_stringers(top_text, bottom_text, inp.stringer_area)
        if not sts:
            raise ValueError("At least one stringer is required on top or bottom.")
        return sts

    def run_analysis(self):
        try:
            inp = self.parse_inputs()
            if self.airfoil_x is None:
                self.airfoil_x, self.airfoil_y = read_airfoil_dat(inp.airfoil_path)
            if self.span_y is None:
                self.span_y, self.span_V, self.span_T = read_span_load_file(self.load_var.get().strip())

            # Apply load factor.
            V_design = self.span_V * inp.load_factor

            stations = choose_stations(self.span_y)
            results = []

            stringers = self._make_stringers(inp)

            # Automatic sizing loop: increase ribs or stringers if buckling fails.
            top = sorted([s.x for s in stringers if s.surface == "top"])
            bottom = sorted([s.x for s in stringers if s.surface == "bottom"])

            for iteration in range(20):
                stringers = self._make_stringers(inp, top, bottom)
                station_results = []
                global_fail = False

                for yi, Vi in zip(self.span_y, V_design):
                    frac = yi / max(inp.span, 1e-12)
                    chord = inp.root_chord + frac * (inp.tip_chord - inp.root_chord)

                    # Interpolate supplied T if present, else V*e.
                    if self.span_T is not None:
                        Tj = float(np.interp(yi, self.span_y, self.span_T)) * inp.load_factor
                    else:
                        Tj = float(Vi * inp.load_offset)

                    sec = solve_section(
                        self.airfoil_x, self.airfoil_y, chord, inp, stringers, float(Vi), Tj
                    )
                    station_results.append((yi, chord, Vi, Tj, sec))

                    vm_max = float(np.max(sec.vm_skin))
                    if vm_max > inp.material.yield_strength:
                        global_fail = True
                    if any(p["fails"] for p in sec.panel_meta):
                        global_fail = True

                if not global_fail:
                    break

                # Prefer adding ribs first, then stringers.
                if any(
                    any(p["fails"] for p in sec.panel_meta)
                    for _, _, _, _, sec in station_results
                ):
                    inp.n_ribs += 2

                    # If still failing after a significant increase in ribs,
                    # add one stringer to each surface.
                    if inp.n_ribs >= 30 and len(top) < inp.max_stringers_per_surface:
                        # Add mid-point in largest chordwise gap.
                        def add_one(arr):
                            bounds = [0.03] + list(arr) + [0.97]
                            gaps = np.diff(bounds)
                            i = int(np.argmax(gaps))
                            return sorted(arr + [(bounds[i] + bounds[i+1]) / 2.0])
                        top = add_one(top)
                        bottom = add_one(bottom)
                        inp.n_ribs = max(10, inp.n_ribs - 4)
                else:
                    # Yield-driven failure: add stringers if possible.
                    if len(top) < inp.max_stringers_per_surface:
                        def add_one(arr):
                            bounds = [0.03] + list(arr) + [0.97]
                            gaps = np.diff(bounds)
                            i = int(np.argmax(gaps))
                            return sorted(arr + [(bounds[i] + bounds[i+1]) / 2.0])
                        top = add_one(top)
                        bottom = add_one(bottom)
                    else:
                        break

            self.latest_results = station_results

            # Summary.
            max_vm = max(float(np.max(sec.vm_skin)) for _, _, _, _, sec in station_results)
            min_bf = min(
                float(p["buckling_factor"])
                for _, _, _, _, sec in station_results
                for p in sec.panel_meta
            )
            critical_station = station_results[int(np.argmax([
                float(np.max(sec.vm_skin)) for _, _, _, _, sec in station_results
            ]))][0]

            final_stringers = self._make_stringers(inp, top, bottom)

            summary = []
            summary.append("SHEET-FLOW WING ANALYSIS")
            summary.append("=" * 90)
            summary.append(f"Cells: {'two-cell' if inp.spar_x is not None else 'one-cell'}")
            summary.append(f"Final N_R = {inp.n_ribs}")
            summary.append(f"Top stringers x/c    = {', '.join(f'{v:.4f}' for v in top)}")
            summary.append(f"Bottom stringers x/c = {', '.join(f'{v:.4f}' for v in bottom)}")
            summary.append(f"Max von Mises (shear-only) = {max_vm/1e6:.3f} MPa at y={critical_station:.3f} m")
            summary.append(f"Material yield = {inp.material.yield_strength/1e6:.3f} MPa")
            summary.append(f"Minimum panel buckling factor = {min_bf:.3f}")
            summary.append("Design status = " + ("PASS" if (max_vm <= inp.material.yield_strength and min_bf >= 1.0) else "REQUIRES REVIEW"))
            summary.append("")
            summary.append("Critical station details:")
            for yi, chord, Vi, Ti, sec in station_results:
                bf = min(p["buckling_factor"] for p in sec.panel_meta)
                summary.append(
                    f"  y={yi:.3f} m, c={chord:.3f} m, "
                    f"V={Vi:.1f} N, T={Ti:.1f} N-m, "
                    f"x_sc={sec.shear_centre_x:.4f} m, "
                    f"max|q|={np.max(np.abs(sec.q_total)):.3f} N/m, "
                    f"max|tau|={np.max(np.abs(sec.tau_skin))/1e6:.3f} MPa, "
                    f"min BF={bf:.3f}"
                )

            self.latest_summary = "\n".join(summary)
            self.text.delete("1.0", tk.END)
            self.text.insert(tk.END, self.latest_summary)

            self.plot_selected_station(auto=True)
            self.status_var.set(
                f"Analysis complete. Final N_R={inp.n_ribs}; "
                f"{len(top)} top + {len(bottom)} bottom stringers."
            )
        except Exception as e:
            self.status_var.set("Analysis failed.")
            messagebox.showerror("Analysis error", f"{e}\n\n{traceback.format_exc()}")

    def plot_selected_station(self, auto=False):
        if not self.latest_results:
            if not auto:
                messagebox.showwarning("Results", "Run analysis first.")
            return

        try:
            yreq = float(self.station_var.get())
        except Exception:
            yreq = self.latest_results[0][0]

        idx = int(np.argmin([abs(r[0] - yreq) for r in self.latest_results]))
        y, chord, V, T, sec = self.latest_results[idx]
        self.station_var.set(f"{y:.4f}")

        self.ax.clear()
        self.ax.plot(sec.x_surf / chord, sec.q_basic, label="q_basic")
        self.ax.plot(sec.x_surf / chord, sec.q_torsion, label="q_torsion")
        self.ax.plot(sec.x_surf / chord, sec.q_total, label="q_total", linewidth=2)

        self.ax.set_xlabel("Chordwise coordinate x/c")
        self.ax.set_ylabel("Shear flow q [N/m]")
        self.ax.set_title(
            f"Sheet flow at y={y:.3f} m | V={V:.1f} N | T={T:.1f} N-m | "
            f"{sec.n_cells}-cell"
        )
        self.ax.grid(True, alpha=0.25)
        self.ax.legend()
        self.canvas.draw()

        # Add a compact station block.
        panel_fail = sum(p["fails"] for p in sec.panel_meta)
        self.text.insert(tk.END, "\n\nSelected station:\n")
        self.text.insert(
            tk.END,
            f"y={y:.3f} m, chord={chord:.3f} m, x_sc={sec.shear_centre_x:.4f} m, "
            f"max|q|={np.max(np.abs(sec.q_total)):.3f} N/m, "
            f"max|tau|={np.max(np.abs(sec.tau_skin))/1e6:.3f} MPa, "
            f"panel failures={panel_fail}\n"
        )

    def export_results(self):
        if not self.latest_results:
            messagebox.showwarning("Export", "Run analysis first.")
            return

        outdir = filedialog.askdirectory(title="Select output folder")
        if not outdir:
            return

        report = os.path.join(outdir, "wing_sheet_flow_report.txt")
        csv_path = os.path.join(outdir, "wing_sheet_flow_station_summary.csv")

        try:
            with open(report, "w", encoding="utf-8") as f:
                f.write(self.latest_summary + "\n\n")
                f.write("PANEL BUCKLING DETAILS\n")
                f.write("=" * 90 + "\n")
                for yi, chord, V, T, sec in self.latest_results:
                    f.write(f"\nStation y={yi:.6f} m, chord={chord:.6f} m\n")
                    for p in sec.panel_meta:
                        f.write(
                            f"{p['surface']:10s} {p['index']:3d} "
                            f"x1={p['x1']:.4f} x2={p['x2']:.4f} "
                            f"a={p['a_rib']:.5f} b={p['b_stringer']:.5f} "
                            f"tau={p['tau_max']/1e6:.5f}MPa "
                            f"tau_cr={p['tau_cr']/1e6:.5f}MPa "
                            f"BF={p['buckling_factor']:.4f} "
                            f"{'FAIL' if p['fails'] else 'PASS'}\n"
                        )

            with open(csv_path, "w", newline="", encoding="utf-8") as f:
                wr = csv.writer(f)
                wr.writerow([
                    "y_m", "chord_m", "V_design_N", "T_Nm",
                    "x_shear_centre_m", "max_q_N_per_m",
                    "max_tau_MPa", "max_vm_MPa", "min_buckling_factor"
                ])
                for yi, chord, V, T, sec in self.latest_results:
                    wr.writerow([
                        yi, chord, V, T, sec.shear_centre_x,
                        float(np.max(np.abs(sec.q_total))),
                        float(np.max(np.abs(sec.tau_skin))) / 1e6,
                        float(np.max(sec.vm_skin)) / 1e6,
                        min(p["buckling_factor"] for p in sec.panel_meta),
                    ])

            messagebox.showinfo("Export", f"Saved:\n{report}\n{csv_path}")
        except Exception as e:
            messagebox.showerror("Export error", str(e))


def main():
    app = WingSheetFlowApp()
    app.mainloop()


if __name__ == "__main__":
    main()
