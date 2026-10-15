"""Simple tip-clearance modelling for small rotors.

Tip gap in micro turbomachinery is set by what can be built and kept, not by what
aerodynamics would like: at 0.15 mm on a 6 mm exit blade the clearance ratio is already 2.5 %.
This module provides

* :func:`disk_rim_growth`, :func:`thermal_growth` - first-order radial growth of a rotor rim
* :func:`clearance_budget` - stack-up of manufacturing, runout, growth and thermal terms
* :func:`clearance_effect` - first-order efficiency and stall-margin penalty for a given gap
* :func:`gap_adjusted_height`, :func:`apply_tip_gap_geometry` - hooks into the geometry builders

All penalty coefficients are generic literature-range defaults. They give the right *trend*
(roughly linear in gap/height) and are meant to be calibrated against rig data.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Literal, Optional
import math


# material data: density kg/m3, Young's modulus Pa, Poisson ratio, CTE 1/K
MATERIALS = {
    "aluminium_7075": dict(rho=2810.0, E=71.7e9, nu=0.33, cte=23.4e-6),
    "titanium_6al4v": dict(rho=4430.0, E=113.8e9, nu=0.342, cte=8.6e-6),
    "inconel_718": dict(rho=8190.0, E=205e9, nu=0.29, cte=13.0e-6),
    "steel_17_4ph": dict(rho=7800.0, E=196e9, nu=0.27, cte=10.8e-6),
}


def disk_rim_growth(radius: float, tip_speed: float, material: str = "titanium_6al4v") -> float:
    """Radial rim growth (m) of a solid rotating disk, ``u = (1 - nu) rho U^2 r / (4 E)``.

    Solid-disk closed form. A bored disk grows somewhat more; blade pull adds a little.
    Temperature effects on E are ignored.
    """
    m = MATERIALS[material]
    return (1.0 - m["nu"]) * m["rho"] * tip_speed**2 * radius / (4.0 * m["E"])


def thermal_growth(radius: float, delta_T: float, material: str = "titanium_6al4v") -> float:
    """Radial thermal growth (m) of a rim at ``delta_T`` above assembly temperature."""
    return MATERIALS[material]["cte"] * delta_T * radius


@dataclass
class ClearanceBudget:
    """Cold-build clearance stack-up. All values in metres."""
    manufacturing: float
    runout: float
    disk_growth: float
    thermal_mismatch: float
    out_of_round: float
    worst_case: float = 0.0
    rss: float = 0.0
    recommended: float = 0.0
    notes: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def clearance_budget(
    radius: float,
    tip_speed: float,
    *,
    material: str = "titanium_6al4v",
    manufacturing: float = 0.05e-3,
    runout: float = 0.03e-3,
    out_of_round: float = 0.03e-3,
    rotor_delta_T: float = 60.0,
    casing_delta_T: float = 20.0,
    casing_material: Optional[str] = None,
    margin: float = 0.05e-3,
) -> ClearanceBudget:
    """Minimum practical *cold* clearance for a rotor of given tip radius and tip speed.

    The rotor grows outward (disk + thermal); the casing grows too, so only the
    *differential* thermal growth enters. The recommended gap is the RSS stack plus a margin,
    never below the worst-case rotor growth that closes the gap in operation.

    Args:
        manufacturing: combined rotor and casing machining tolerance on the gap.
        runout: shaft/bearing radial runout (dynamic).
        out_of_round: casing and wheel out-of-roundness.
        rotor_delta_T, casing_delta_T: temperature rise of rotor rim and casing above build.
    """
    grow_disk = disk_rim_growth(radius, tip_speed, material)
    grow_rotor_th = thermal_growth(radius, rotor_delta_T, material)
    grow_case_th = thermal_growth(radius, casing_delta_T, casing_material or material)
    thermal_mismatch = max(grow_rotor_th - grow_case_th, 0.0)
    terms = [manufacturing, runout, out_of_round, grow_disk, thermal_mismatch]
    worst = sum(terms)
    rss = math.sqrt(sum(t * t for t in terms))
    # operating closure the gap must at least absorb (deterministic growth terms)
    closure = grow_disk + thermal_mismatch
    recommended = max(rss + margin, closure + margin)
    notes = []
    if grow_disk > 0.5 * recommended:
        notes.append("Disk growth dominates the stack: lower tip speed or use a stiffer material.")
    if thermal_mismatch > 0.5 * recommended:
        notes.append("Thermal mismatch dominates: improve casing/rotor temperature matching.")
    return ClearanceBudget(manufacturing, runout, grow_disk, thermal_mismatch, out_of_round,
                           worst_case=worst, rss=rss, recommended=recommended, notes=notes)


@dataclass
class ClearanceEffect:
    gap: float
    height: float
    ratio_percent: float
    delta_efficiency_points: float
    delta_stall_margin_points: float
    status: str
    note: str


# generic first-order coefficients: efficiency points lost per 1 % of gap/height,
# and stall-margin points lost per 1 % of gap/height
_COEFFS = {
    "axial_compressor": (2.0, 5.0),
    "axial_turbine": (1.8, 0.0),
    "centrifugal": (1.5, 4.0),
}


def clearance_effect(
    gap: float,
    height: float,
    kind: Literal["axial_compressor", "axial_turbine", "centrifugal"] = "centrifugal",
    k_eta: Optional[float] = None,
    k_stall: Optional[float] = None,
) -> ClearanceEffect:
    """First-order effect of a tip gap on efficiency and stall margin.

    ``height`` is the blade height (axial rotors) or exit blade height b2 (centrifugal).
    Efficiency drops by ``k_eta`` points per 1 % of gap/height; defaults are typical literature
    magnitudes, not measured values.
    """
    ke, ks = _COEFFS[kind]
    ke = ke if k_eta is None else k_eta
    ks = ks if k_stall is None else k_stall
    ratio = 100.0 * gap / height
    status = "ok" if ratio <= 2.0 else ("warn" if ratio <= 4.0 else "high")
    note = {
        "ok": "Clearance ratio within a normal range.",
        "warn": "Clearance ratio is high; efficiency and surge margin will be noticeably reduced.",
        "high": "Clearance ratio is very high for this height; consider a larger rotor or a shrouded design.",
    }[status]
    return ClearanceEffect(gap, height, ratio, ke * ratio, ks * ratio, status, note)


def gap_adjusted_height(passage_height: float, gap: float) -> float:
    """Physical blade height for an axial row in a passage of ``passage_height``: ``h - gap``."""
    if gap >= passage_height:
        raise ValueError("gap must be smaller than the passage height")
    return passage_height - gap


def apply_tip_gap_geometry(centrif, gap: float) -> float:
    """Pull the tip of a *built* :class:`microbladegen.aero.Centrif` blade ``gap`` away from the shroud.

    ``Centrif.tip_clearance`` only reparameterises the span variable: the tip profile is still
    evaluated on the shroud, so blade coordinates do not move. This function does the geometric
    part: at every chord station each spanwise column is contracted toward the hub along the
    hub-to-tip line by ``(L - gap) / L`` (``L`` = local hub-to-shroud distance), so the tip ends
    exactly ``gap`` below the casing (measured along the hub-to-shroud cut line, so the
    perpendicular distance to a curved shroud is slightly smaller) and interior sections keep
    their relative positions.
    Cylindrical and Cartesian point arrays of the main blade and splitter are updated in place;
    meridional (mp, theta) arrays are not.

    Returns the largest relative contraction applied (``gap / L_min``).
    """
    import numpy as np
    if gap <= 0:
        return 0.0
    worst = 0.0
    for blade in (centrif.mainblade, getattr(centrif, "splitterblade", None)):
        if blade is None:
            continue
        for cyl, cart in ((blade.ss_cyl_pts, blade.ss_cart_pts), (blade.ps_cyl_pts, blade.ps_cart_pts)):
            hub = cyl[0, :, :2].copy()                          # (nchord, 2) -> x, r
            tip = cyl[-1, :, :2].copy()
            length = np.hypot(*(tip - hub).T)
            f = np.clip((length - gap) / np.maximum(length, 1e-12), 0.2, 1.0)
            worst = max(worst, float(np.max(1.0 - f)))
            cyl[:, :, :2] = hub[None, :, :] + f[None, :, None] * (cyl[:, :, :2] - hub[None, :, :])
            r, th = cyl[:, :, 1], cyl[:, :, 2]
            cart[:, :, 0] = cyl[:, :, 0]
            cart[:, :, 1] = r * np.sin(th)
            cart[:, :, 2] = r * np.cos(th)
    return worst
