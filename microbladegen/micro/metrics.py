"""Blade-row and impeller metrics for micro-turbomachinery.

Axial rows (:func:`axial_row_metrics`) work on a built
:class:`microbladegen.aero.Airfoil3D` (x axial, y tangential, z radial offset from the hub).
Impellers (:func:`impeller_metrics`) work on a :class:`~microbladegen.micro.impeller.MicroImpeller`.

Reported quantities include solidity, aspect ratio, thickness-to-chord, blade speed, inlet and
exit Mach numbers in the blade frame, tip blade Mach number, chord Reynolds number and, for
decelerating rows, the de Haller number.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from microbladegen.micro.gas import Gas, AIR


@dataclass
class RowFlow:
    """Mean-line flow through one blade row (angles in degrees from axial, in the blade frame).

    Attributes:
        rpm: shaft speed (0 for a stator).
        hub_radius: radius of the z = 0 surface (m).
        nblades: blades in the row.
        axial_velocity: Cx (m/s), constant through the row.
        inlet_angle: absolute flow angle at the *mean* radius entering the row. Swirl follows a
            free vortex (``Ctheta r = const``) across the span.
        exit_angle: absolute flow angle at mean radius leaving the row (None to skip exit metrics).
        Tt_in, Pt_in: inlet total temperature (K) and pressure (Pa) in the absolute frame.
        Tt_out: exit total temperature (K); defaults to ``Tt_in`` (stators, turbine NGVs).
        gas: gas model (use ``Gas(gamma=1.33)`` for hot-section air/combustion products).
    """
    rpm: float
    hub_radius: float
    nblades: int
    axial_velocity: float
    inlet_angle: float
    Tt_in: float = 288.15
    Pt_in: float = 101325.0
    exit_angle: Optional[float] = None
    Tt_out: Optional[float] = None
    gas: Gas = AIR


def _max_thickness(ss: np.ndarray, ps: np.ndarray) -> float:
    """Maximum section thickness: largest nearest-neighbour distance from SS to the PS polyline."""
    p, q = ss[:, :2], ps[:, :2]
    d = np.min(np.linalg.norm(p[:, None, :] - q[None, :, :], axis=2), axis=1)
    return float(d.max())


def axial_row_metrics(blade, flow: RowFlow, span_fractions: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0)) -> pd.DataFrame:
    """Per-section metrics of a built :class:`Airfoil3D` row (units: mm, m/s, dimensionless)."""
    g = flow.gas
    chord, _ = blade.get_chord()
    zz = np.asarray(blade.zz, dtype=float)
    span = float(zz[-1] - zz[0])
    r_mean = flow.hub_radius + 0.5 * span
    omega = 2.0 * np.pi * flow.rpm / 60.0
    rotating = flow.rpm > 0
    cx = flow.axial_velocity
    ct1_mean = cx * np.tan(np.radians(flow.inlet_angle))
    ct2_mean = cx * np.tan(np.radians(flow.exit_angle)) if flow.exit_angle is not None else None
    Tt_out = flow.Tt_out if flow.Tt_out is not None else flow.Tt_in
    chord_mean = float(np.mean(chord))

    rows = []
    for f in span_fractions:
        j = int(round(f * (len(zz) - 1)))
        r = flow.hub_radius + zz[j]
        U = omega * r
        ct1 = ct1_mean * r_mean / r                               # free vortex
        c1 = np.hypot(cx, ct1)
        w1 = np.hypot(cx, U - ct1) if rotating else c1
        T1, p1, rho1 = g.static_from_total(flow.Tt_in, flow.Pt_in, c1)
        a1 = float(g.speed_of_sound(T1))
        mu1 = float(g.viscosity(T1))
        pitch = 2.0 * np.pi * r / flow.nblades
        tmax = _max_thickness(blade.shft_ss[j], blade.shft_ps[j])
        row = {
            "span": float(f), "radius_mm": r * 1e3, "chord_mm": float(chord[j]) * 1e3, "pitch_mm": pitch * 1e3,
            "solidity": float(chord[j] / pitch), "aspect_ratio": span / chord_mean,
            "t_max_mm": tmax * 1e3, "t_over_c": tmax / float(chord[j]),
            "blade_speed": U, "inlet_velocity": w1,
            "inlet_mach": w1 / a1, "tip_blade_mach": U / a1 if rotating else 0.0,
            "reynolds_chord": float(rho1 * w1 * chord[j] / mu1),
        }
        if ct2_mean is not None:
            ct2 = ct2_mean * r_mean / r
            c2 = np.hypot(cx, ct2)
            w2 = np.hypot(cx, U - ct2) if rotating else c2
            T2, _, rho2 = g.static_from_total(Tt_out, flow.Pt_in, c2)
            a2 = float(g.speed_of_sound(T2))
            ratio = w2 / w1
            row["exit_velocity"] = w2
            row["exit_mach"] = w2 / a2
            row["reynolds_chord_exit"] = float(rho2 * w2 * chord[j] / float(g.viscosity(T2)))
            row["velocity_ratio"] = ratio
            row["de_haller"] = ratio if ratio < 1.0 else float("nan")     # diffusing rows only
            row["turning_deg"] = float(abs(np.degrees(np.arctan2((U - ct2) if rotating else ct2, cx)
                                                      - np.arctan2((U - ct1) if rotating else ct1, cx))))
        rows.append(row)
    return pd.DataFrame(rows).set_index("span")


def free_vortex_angles(rpm: float, radii: Sequence[float], mean_radius: float, axial_velocity: float,
                       alpha1_mean: float, alpha2_mean: float) -> pd.DataFrame:
    """Flow angles (deg from axial, magnitudes) across the span for a free-vortex row.

    ``alpha1_mean`` / ``alpha2_mean`` are the absolute flow angles entering and leaving the row at
    ``mean_radius``. Returns absolute (alpha) and relative (beta) angles at each radius. For a
    stator (``rpm = 0``) the beta columns equal the alpha columns. Use the result as the blade
    angles (flow angle = metal angle; incidence and deviation are not applied).
    """
    r = np.asarray(radii, dtype=float)
    omega = 2.0 * np.pi * rpm / 60.0
    u = omega * r
    ta1 = np.tan(np.radians(alpha1_mean)) * mean_radius / r
    ta2 = np.tan(np.radians(alpha2_mean)) * mean_radius / r
    if rpm > 0:
        tb1 = np.abs(u / axial_velocity - ta1)
        tb2 = np.abs(u / axial_velocity - ta2)
    else:
        tb1, tb2 = np.abs(ta1), np.abs(ta2)
    return pd.DataFrame({
        "radius_mm": r * 1e3, "U": u,
        "alpha1": np.degrees(np.arctan(np.abs(ta1))), "alpha2": np.degrees(np.arctan(np.abs(ta2))),
        "beta1": np.degrees(np.arctan(tb1)), "beta2": np.degrees(np.arctan(tb2)),
    })


def impeller_metrics(imp) -> pd.DataFrame:
    """Hub / mid / tip geometry and flow metrics of a built impeller (mean-line sizing + geometry)."""
    s = imp.sizing
    mb = imp.centrif.mainblade
    n = mb.ss_cyl_pts.shape[0]
    rows = []
    for name, j, f in (("hub", 0, 0.0), ("mid", n // 2, 0.5), ("tip", n - 1, 1.0)):
        cam = 0.5 * (mb.ss_cart_pts[j] + mb.ps_cart_pts[j])
        length = float(np.sum(np.linalg.norm(np.diff(cam, axis=0), axis=1)))
        thick = float(np.max(np.linalg.norm(mb.ss_cart_pts[j] - mb.ps_cart_pts[j], axis=1)))
        pitch_exit_full = 2.0 * np.pi * s.r2 / s.n_full
        rows.append({
            "section": name, "span": f, "blade_length_mm": length * 1e3, "t_max_mm": thick * 1e3,
            "t_over_length": thick / length,
            "solidity_exit_full_blades": length / pitch_exit_full,
            "solidity_exit_incl_splitters": length * (s.n_full + 0.5 * s.n_splitter) / (2.0 * np.pi * s.r2),
            "inlet_angle_deg": imp.measured["inlet_angle_hub_mid_tip"][("hub", "mid", "tip").index(name)],
            "exit_angle_deg": imp.measured["exit_angle_hub_mid_tip"][("hub", "mid", "tip").index(name)],
            "wrap_deg": imp.measured["wrap_angle_hub_mid_tip_deg"][("hub", "mid", "tip").index(name)],
        })
    return pd.DataFrame(rows).set_index("section")


def impeller_flow_table(s) -> pd.Series:
    """One-column table of the mean-line impeller quantities a micro-turbojet designer reviews first."""
    d = {
        "Mass flow (kg/s)": s.mass_flow, "Speed (rpm)": s.rpm, "Impeller PR (tt)": s.pressure_ratio,
        "Efficiency (tt)": s.eta_tt, "Tip radius r2 (mm)": s.r2 * 1e3, "Exit blade height b2 (mm)": s.b2 * 1e3,
        "Tip speed U2 (m/s)": s.u2, "Tip Mach U2/a01": s.tip_mach, "Machine Reynolds rho01 U2 r2 / mu01": s.machine_reynolds,
        "Slip factor (Wiesner)": s.slip_factor, "Work coefficient dh0/U2^2": s.work_coefficient,
        "Exit flow coefficient cr2/U2": s.cr2 / s.u2, "Absolute exit Mach": s.exit_mach_abs,
        "Backsweep from radial (deg)": s.backsweep, "Exit relative flow angle (deg)": s.beta2_flow,
        "Inducer shroud radius r1s (mm)": s.r1s * 1e3, "Inducer hub radius r1h (mm)": s.r1h * 1e3,
        "r1s / r2": s.radius_ratio, "Shroud relative Mach W1s/a1": s.inlet_rel_mach_shroud,
        "Shroud flow angle (deg)": s.beta1s_flow, "Diffusion ratio W1s/W2": s.diffusion_ratio,
        "Blades (full + splitter)": f"{s.n_full} + {s.n_splitter}", "Specific speed Ns": s.specific_speed,
        "Axial length (mm)": s.axial_length * 1e3,
    }
    return pd.Series(d, name="value")
