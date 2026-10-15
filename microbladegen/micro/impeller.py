"""Centrifugal impeller design for micro-turbojets.

Two layers:

* :func:`size_impeller` - mean-line sizing from engine requirements
  (mass flow, speed, impeller pressure ratio). Returns tip speed, tip Mach,
  inducer and exit dimensions, blade count, diffusion ratio, Reynolds number.
* :func:`build_impeller` - turns the sizing into a 3D blade through the
  existing :class:`microbladegen.aero.Centrif` generator. The generator's
  nominal metal-angle inputs do not map one-to-one to the angles that come
  out, so the builder *calibrates* wrap and leading-edge inputs until the
  measured inlet and exit blade angles match the sizing targets.

Conventions: angles in degrees. Blade inlet angle is measured from the axial
(meridional) direction; exit backsweep from the radial direction.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np

from microbladegen.micro.gas import Gas, AIR
from microbladegen.micro.tip_clearance import apply_tip_gap_geometry


# --------------------------------------------------------------------------- #
# Sizing
# --------------------------------------------------------------------------- #
@dataclass
class ImpellerSizing:
    """Mean-line result. Lengths in m, speeds in m/s, angles in deg."""
    # requirements
    mass_flow: float
    rpm: float
    pressure_ratio: float
    eta_tt: float
    Tt1: float
    Pt1: float
    # exit
    r2: float
    b2: float
    u2: float
    tip_mach: float                 # U2 / a01
    machine_reynolds: float         # rho01 U2 r2 / mu01
    slip_factor: float
    work_coefficient: float         # dh0 / U2^2
    specific_work: float            # J/kg
    cr2: float
    ctheta2: float
    w2: float
    beta2_flow: float               # relative flow angle at exit, from radial
    backsweep: float                # blade exit angle, from radial
    exit_mach_abs: float
    # inlet
    r1s: float
    r1h: float
    hub_to_shroud: float
    cx1: float
    beta1s_flow: float              # shroud relative flow angle, from axial
    beta1s_blade: float             # shroud blade angle (flow minus incidence)
    beta1h_blade: float
    w1s: float
    inlet_rel_mach_shroud: float
    # blading
    n_full: int
    n_splitter: int
    diffusion_ratio: float          # W1s / W2
    radius_ratio: float             # r1s / r2
    specific_speed: float
    axial_length: float
    notes: list = field(default_factory=list)

    @property
    def n_total(self) -> int:
        return self.n_full + self.n_splitter

    def as_dict(self) -> dict:
        return asdict(self)


def wiesner_slip(n_blades: int, backsweep_deg: float) -> float:
    """Wiesner slip factor, sigma = 1 - sqrt(cos(beta2b)) / Z^0.7 (beta2b from radial)."""
    return 1.0 - np.sqrt(np.cos(np.radians(backsweep_deg))) / n_blades**0.7


def glassman_blade_count(alpha2_deg: float) -> float:
    """Glassman (1976) blade-count heuristic, Z = (pi/30)(110 - alpha2) tan(alpha2).

    alpha2 is the absolute exit flow angle from radial, in degrees.
    """
    a = np.radians(alpha2_deg)
    return float(np.pi / 30.0 * (110.0 - alpha2_deg) * np.tan(a))


def size_impeller(
    mass_flow: float,
    rpm: float,
    pressure_ratio: float,
    *,
    eta_tt: float = 0.82,
    Tt1: float = 288.15,
    Pt1: float = 101325.0,
    gas: Gas = AIR,
    backsweep: float = 30.0,
    beta1s_flow: float = 58.0,
    incidence: float = 3.0,
    hub_to_shroud: float = 0.36,
    exit_flow_coefficient: float = 0.25,
    disk_friction: float = 0.05,
    inlet_blockage: float = 0.05,
    exit_blockage: float = 0.10,
    n_full: Optional[int] = None,
    splitters: bool = True,
    axial_length_ratio: float = 0.45,
) -> ImpellerSizing:
    """Mean-line sizing of a micro-turbojet centrifugal impeller.

    Args:
        mass_flow: impeller mass flow (kg/s).
        rpm: shaft speed.
        pressure_ratio: total-to-total pressure ratio across the *impeller* (exit of rotor).
        eta_tt: impeller total-to-total efficiency.
        backsweep: blade exit angle from radial (deg), 0 = radial.
        beta1s_flow: target relative flow angle at the inducer shroud, from axial (deg).
            50-62 deg keeps shroud relative Mach near its minimum for a given flow.
        incidence: blade angle = flow angle - incidence (deg).
        hub_to_shroud: inducer r1h / r1s.
        exit_flow_coefficient: cr2/U2.
        disk_friction: windage/disk-friction work as a fraction of Euler work.
        n_full, splitters: full-blade count; if None a Glassman estimate is used (rounded to
            an even total). With ``splitters=True`` the total is split evenly into full blades
            and splitter blades.
        axial_length_ratio: axial length / r2 used for the meridional contour.
    """
    g = gas
    cp, gam, R = g.cp, g.gamma, g.R
    w = 2.0 * np.pi * rpm / 60.0
    notes: list[str] = []

    dh_is = cp * Tt1 * (pressure_ratio ** ((gam - 1.0) / gam) - 1.0)
    dh0 = dh_is / eta_tt

    # blade count and slip are coupled; iterate a few times
    n_total = 12
    sigma = wiesner_slip(n_total, backsweep)
    for _ in range(12):
        mu = sigma - exit_flow_coefficient * np.tan(np.radians(backsweep))   # ctheta2/U2
        u2 = np.sqrt(dh0 / ((1.0 + disk_friction) * mu))
        alpha2 = np.degrees(np.arctan2(mu, exit_flow_coefficient))
        if n_full is None:
            n_est = glassman_blade_count(alpha2)
            n_new = max(8, int(2 * round(n_est / 2)))
        else:
            n_new = 2 * n_full if splitters else n_full
        if n_new == n_total:
            break
        n_total = n_new
        sigma = wiesner_slip(n_total, backsweep)
    n_f = n_full if n_full is not None else (n_total // 2 if splitters else n_total)
    n_s = n_f if splitters else 0
    n_total = n_f + n_s
    sigma = wiesner_slip(n_total, backsweep)
    mu = sigma - exit_flow_coefficient * np.tan(np.radians(backsweep))
    u2 = np.sqrt(dh0 / ((1.0 + disk_friction) * mu))
    r2 = u2 / w
    cr2 = exit_flow_coefficient * u2
    ct2 = mu * u2

    # inlet stagnation reference
    a01 = float(g.speed_of_sound(Tt1))
    rho01 = Pt1 / (R * Tt1)
    mu01 = float(g.viscosity(Tt1))
    mach_u2 = u2 / a01
    re_m = rho01 * u2 * r2 / mu01

    # exit static state (total pressure from PR; total temperature from work)
    Tt2 = Tt1 + dh0 / cp
    Pt2 = Pt1 * pressure_ratio
    c2 = np.hypot(cr2, ct2)
    T2, p2, rho2 = g.static_from_total(Tt2, Pt2, c2)
    b2 = mass_flow / (rho2 * cr2 * 2.0 * np.pi * r2 * (1.0 - exit_blockage))
    wt2 = u2 - ct2
    w2 = np.hypot(cr2, wt2)
    beta2_flow = float(np.degrees(np.arctan2(wt2, cr2)))
    a2 = float(g.speed_of_sound(T2))

    # inducer: solve r1s from continuity at the chosen shroud flow angle
    tb = np.tan(np.radians(beta1s_flow))
    nu = hub_to_shroud
    r1s = 0.5 * r2
    for _ in range(100):
        u1s = w * r1s
        cx1 = u1s / tb
        T1, p1, rho1 = g.static_from_total(Tt1, Pt1, cx1)
        area = mass_flow / (rho1 * cx1 * (1.0 - inlet_blockage))
        r_new = np.sqrt(area / (np.pi * (1.0 - nu**2)))
        if abs(r_new - r1s) < 1e-10:
            r1s = r_new
            break
        r1s = 0.5 * r1s + 0.5 * r_new
    u1s = w * r1s
    cx1 = u1s / tb
    T1, p1, rho1 = g.static_from_total(Tt1, Pt1, cx1)
    w1s = np.hypot(cx1, u1s)
    a1 = float(g.speed_of_sound(T1))
    m_w1s = w1s / a1
    r1h = nu * r1s
    beta1h_flow = np.degrees(np.arctan(nu * tb))     # constant cx1 across span

    # diffusion ratio and specific speed
    dr = w1s / w2
    q1 = mass_flow / rho01
    ns = w * np.sqrt(q1) / dh_is**0.75

    # advisory notes
    if mach_u2 > 1.35:
        notes.append(f"Tip Mach {mach_u2:.2f} is high for an unshrouded micro impeller; expect supersonic diffuser entry.")
    if m_w1s > 1.2:
        notes.append(f"Shroud relative Mach {m_w1s:.2f} exceeds ~1.2; consider larger r1s (lower beta1s) or more speed margin.")
    elif m_w1s > 0.95:
        notes.append(f"Shroud relative Mach {m_w1s:.2f}: transonic inducer.")
    if dr > 1.9:
        notes.append(f"Diffusion ratio {dr:.2f} > 1.9: risk of impeller stall; reduce backsweep or increase r2.")
    if b2 < 0.0025:
        notes.append(f"Exit blade height {b2*1e3:.2f} mm is very small: tip clearance will dominate losses.")
    if u2 > 480.0:
        notes.append(f"Tip speed {u2:.0f} m/s: check disk stress for the intended material.")

    return ImpellerSizing(
        mass_flow=mass_flow, rpm=rpm, pressure_ratio=pressure_ratio, eta_tt=eta_tt, Tt1=Tt1, Pt1=Pt1,
        r2=float(r2), b2=float(b2), u2=float(u2), tip_mach=float(mach_u2), machine_reynolds=float(re_m),
        slip_factor=float(sigma), work_coefficient=float(dh0 / u2**2), specific_work=float(dh0),
        cr2=float(cr2), ctheta2=float(ct2), w2=float(w2), beta2_flow=beta2_flow, backsweep=float(backsweep),
        exit_mach_abs=float(c2 / a2),
        r1s=float(r1s), r1h=float(r1h), hub_to_shroud=float(nu), cx1=float(cx1), beta1s_flow=float(beta1s_flow),
        beta1s_blade=float(beta1s_flow - incidence), beta1h_blade=float(beta1h_flow - incidence),
        w1s=float(w1s), inlet_rel_mach_shroud=float(m_w1s),
        n_full=int(n_f), n_splitter=int(n_s), diffusion_ratio=float(dr), radius_ratio=float(r1s / r2),
        specific_speed=float(ns), axial_length=float(axial_length_ratio * r2), notes=notes,
    )


# --------------------------------------------------------------------------- #
# Meridional contour and 3D build
# --------------------------------------------------------------------------- #
def meridional_endwalls(s: ImpellerSizing, npts: int = 120):
    """Hub and shroud (x, r) curves: elliptical axial-to-radial turn.

    The shroud ends at x = L - b2 and the hub (back plate) at x = L, both with radial tangent at r2.
    Returns ``(hub, shroud)`` as (npts, 2) arrays.
    """
    th = np.linspace(0.0, np.pi / 2.0, npts)
    L, b2 = s.axial_length, s.b2
    shroud = np.column_stack([(L - b2) * np.sin(th), s.r2 - (s.r2 - s.r1s) * np.cos(th)])
    hub = np.column_stack([L * np.sin(th), s.r2 - (s.r2 - s.r1h) * np.cos(th)])
    return hub, shroud


def _arc_length(curve: np.ndarray) -> float:
    return float(np.sum(np.hypot(np.diff(curve[:, 0]), np.diff(curve[:, 1]))))


@dataclass
class MicroImpeller:
    """A built impeller: sizing, endwalls, the Centrif object and measured geometry."""
    sizing: ImpellerSizing
    centrif: object
    hub: np.ndarray
    shroud: np.ndarray
    tip_clearance: float
    wraps: list
    le_inputs: list
    measured: dict

    @property
    def nblades(self) -> int:
        return self.sizing.n_full


def _mean_camber_angles(blade, j: int):
    c = 0.5 * (blade.ss_cyl_pts[j] + blade.ps_cyl_pts[j])
    x, r, t = c[:, 0], c[:, 1], c[:, 2]
    dm = np.hypot(np.diff(x), np.diff(r))
    ang = np.degrees(np.arctan2(0.5 * (r[1:] + r[:-1]) * np.diff(t), dm))
    return ang, t


def build_impeller(
    sizing: ImpellerSizing,
    *,
    tip_clearance: float = 0.0,
    tip_thickness: float = 0.7e-3,
    root_thickness_factor: float = 1.5,
    npts_span: int = 11,
    npts_chord: int = 100,
    wrap_bias: float = 0.0,
    inlet_station: int = 8,
    exit_station: int = -8,
    refine_passes: int = 4,
) -> MicroImpeller:
    """Build a 3D impeller (main blades + splitters) from a sizing result.

    The blade camber is a cubic Bezier in the (m', theta) plane, ``m' = integral(dm / r)``, where
    the curve slope is exactly ``tan(beta)``. The builder therefore uses the non-ray camber mode,
    in which the end tangents equal ``LE_Metal_Angle`` and ``TE_Metal_Angle``. Both are given with
    the same sign so the blade leans one way from inducer to exit (no S-curve).

    Args:
        tip_clearance: gap at the shroud (m). Applied geometrically after the build by
            :func:`~microbladegen.micro.tip_clearance.apply_tip_gap_geometry`.
        tip_thickness: blade thickness at the tip section (m). ``Centrif`` takes thickness as a
            fraction of camber length per side, so the builder converts; the root is thicker by
            ``root_thickness_factor`` and the mid span in between.
        wrap_bias: degrees added to the mean-angle estimate of the wrap input (shifts loading
            toward the inducer when positive, toward the exit when negative).
        inlet_station, exit_station: chord indices where the measured angles are reported.
        refine_passes: extra builds that nudge the LE/TE angle inputs by the measured error
            (each build takes about a second).
    """
    from microbladegen.aero import Centrif, CentrifProfile, TrailingEdgeProperties

    hub, shroud = meridional_endwalls(sizing)
    length = 0.5 * (_arc_length(hub) + _arc_length(shroud))
    f_tip = tip_thickness / (2.0 * length)
    shape = np.array([1.0, 1.0, 0.9, 0.75, 0.6])
    spans = (0.0, 0.5, 1.0)

    tb_s = np.tan(np.radians(sizing.beta1s_flow))
    inc = sizing.beta1s_flow - sizing.beta1s_blade
    nu_mid = 0.5 * (1.0 + sizing.hub_to_shroud)
    flow_angles = [np.degrees(np.arctan(sizing.hub_to_shroud * tb_s)),
                   np.degrees(np.arctan(nu_mid * tb_s)), sizing.beta1s_flow]
    le_targets = [fa - inc for fa in flow_angles]            # blade angle from axial, magnitude
    te_target = abs(sizing.backsweep)
    thick = [f_tip * root_thickness_factor, f_tip * (1.0 + 0.5 * (root_thickness_factor - 1.0)), f_tip]

    def make(le_in, te_in):
        cen = Centrif(blade_position=(0.0, 1.0), use_mid_wrap_angle=False, use_ray_camber=False)
        cen.add_hub(hub[:, 0], hub[:, 1])
        cen.add_shroud(shroud[:, 0], shroud[:, 1])
        te = TrailingEdgeProperties(TE_Cut=False, TE_Radius=0.05)
        profs, wr = [], []
        for sp, le, te_a, tk in zip(spans, le_in, te_in, thick):
            wrap = -(0.5 * (le + te_a) + wrap_bias)
            wr.append(wrap)
            profs.append(CentrifProfile(
                percent_span=sp, LE_Thickness=tk, trailing_edge_properties=te,
                LE_Metal_Angle=-le, TE_Metal_Angle=-te_a,
                LE_Metal_Angle_Loc=0.1, TE_Metal_Angle_Loc=0.9,
                ss_thickness=list(tk * shape), ps_thickness=list(tk * shape),
                wrap_angle=wrap, use_ray_camber=False))
        for p in profs:
            cen.add_profile(p)
        if sizing.n_splitter:
            cen.add_splitter(profs, splitter_start=0.5)
        cen.tip_clearance = tip_clearance     # keeps upstream streamline helpers consistent
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cen.build(npts_span=npts_span, npts_chord=npts_chord, nblades=sizing.n_full,
                      nsplitters=1 if sizing.n_splitter else 0)
        # upstream only reparameterises span; move the blade geometry itself
        apply_tip_gap_geometry(cen, tip_clearance)
        return cen, wr

    def measure(cen):
        mb_ = cen.mainblade
        n = mb_.ss_cyl_pts.shape[0]
        out = []
        for j in (0, n // 2, n - 1):
            ang, t = _mean_camber_angles(mb_, j)
            out.append((float(abs(ang[inlet_station])), float(abs(ang[exit_station])), float(np.degrees(abs(t[-1] - t[0])))))
        return out

    le_in = list(le_targets)
    te_in = [te_target] * 3
    cen, wraps = make(le_in, te_in)
    for _ in range(refine_passes):
        m = measure(cen)
        le_in = [float(np.clip(li + 0.8 * (lt - mm[0]), 5.0, 80.0)) for li, lt, mm in zip(le_in, le_targets, m)]
        te_in = [float(np.clip(ti + 0.8 * (te_target - mm[1]), 5.0, 80.0)) for ti, mm in zip(te_in, m)]
        cen, wraps = make(le_in, te_in)

    mb = cen.mainblade
    js = (0, mb.ss_cyl_pts.shape[0] // 2, mb.ss_cyl_pts.shape[0] - 1)
    m = measure(cen)
    le_m = [v[0] for v in m]
    te_m = [v[1] for v in m]
    wrap_m = [v[2] for v in m]
    thick_m = [float(np.max(np.linalg.norm(mb.ss_cart_pts[j] - mb.ps_cart_pts[j], axis=1))) for j in js]
    measured = {
        "inlet_angle_target_hub_mid_tip": [float(v) for v in le_targets],
        "inlet_angle_hub_mid_tip": le_m,
        "exit_angle_target": float(te_target),
        "exit_angle_hub_mid_tip": te_m,
        "max_thickness_hub_mid_tip": thick_m,
        "wrap_angle_hub_mid_tip_deg": wrap_m,
        "exit_radius": float(np.max(mb.ss_cyl_pts[:, :, 1])),
        "camber_length": float(length),
    }
    return MicroImpeller(sizing=sizing, centrif=cen, hub=hub, shroud=shroud, tip_clearance=tip_clearance,
                         wraps=wraps, le_inputs=[-v for v in le_in], measured=measured)
