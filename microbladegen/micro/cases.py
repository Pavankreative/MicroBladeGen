"""Ready-made micro-turbojet design cases.

A case is a JSON file (see ``microbladegen/micro/cases/``) with the requirements or section
definitions of one blade row or impeller. :func:`build_case` turns it into geometry plus
metrics; the same JSON can be copied and edited to make your own case.

Case types
----------
``centrifugal_impeller``
    Mean-line requirements -> :func:`size_impeller` -> :func:`build_impeller`.
``axial_row``
    Hub/mid/tip sections from free-vortex triangles -> :class:`Airfoil3D` with stacking, lean and
    sweep -> :func:`axial_row_metrics`.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib import resources
from typing import Any, Optional

import numpy as np
import pandas as pd

from microbladegen.micro.gas import Gas
from microbladegen.micro.impeller import ImpellerSizing, MicroImpeller, size_impeller, build_impeller
from microbladegen.micro.metrics import RowFlow, axial_row_metrics, free_vortex_angles, impeller_metrics
from microbladegen.micro.tip_clearance import clearance_budget, clearance_effect, ClearanceBudget, ClearanceEffect

_PKG = "microbladegen.micro"


def list_cases() -> list[str]:
    """Names of the bundled cases."""
    folder = resources.files(_PKG).joinpath("cases")
    return sorted(p.name[:-5] for p in folder.iterdir() if p.name.endswith(".json"))


def load_case(name_or_path: str) -> dict:
    """Load a bundled case by name, or any case JSON by path."""
    if name_or_path.endswith(".json"):
        with open(name_or_path) as f:
            return json.load(f)
    path = resources.files(_PKG).joinpath("cases", f"{name_or_path}.json")
    if not path.is_file():
        raise FileNotFoundError(f"unknown case '{name_or_path}'. Available: {', '.join(list_cases())}")
    return json.loads(path.read_text())


@dataclass
class ImpellerCase:
    case: dict
    sizing: ImpellerSizing
    impeller: MicroImpeller
    clearance_budget: ClearanceBudget
    clearance_effect: ClearanceEffect
    metrics: pd.DataFrame


@dataclass
class AxialRowCase:
    case: dict
    blade: Any                       # Airfoil3D
    flow: RowFlow
    angles: pd.DataFrame
    metrics: pd.DataFrame
    nblades: int
    hub_radius: float
    tip_radius: float
    tip_gap: float = 0.0


def _make_foil(alpha1, alpha2, axial_chord, stagger, th, pitch_to_chord=0.75):
    from microbladegen.aero import Airfoil2D
    f = Airfoil2D(alpha1=float(alpha1), alpha2=float(alpha2), axial_chord=float(axial_chord), stagger=float(stagger))
    f.add_le_thickness(th["le"])
    f.add_ps_thickness(thicknessArray=list(th["ps"]), expansion_ratio=1.2)
    f.add_ss_thickness(thicknessArray=list(th["ss"]), camberPercent=0.8, expansion_ratio=1.2)
    f.match_le_thickness()
    f.te_create(radius=th["te_radius_mm"] * 1e-3, wedge_ss=3.0, wedge_ps=2.4)
    f.add_ss_flow_guidance_2(s_c=float(pitch_to_chord), n=10)   # s/c sets the throat the SS straightening aims at
    return f


def axial_row_foils(case: dict):
    """Hub/mid/tip :class:`Airfoil2D` sections of an ``axial_row`` case.

    Returns ``(foils, spans, angles, height, hub_radius, tip_radius, tip_gap)`` where ``height``
    is the physical blade height (passage height minus tip gap) and ``angles`` the free-vortex
    flow-angle table the sections were built from.
    """
    row, sec, th = case["row"], case["sections"], case["thickness"]
    rh, rt = row["hub_radius_mm"] * 1e-3, row["tip_radius_mm"] * 1e-3
    tip_gap = case.get("tip_gap_mm", 0.0) * 1e-3
    height = rt - rh - tip_gap
    rpm = row["rpm"] if row.get("rotating", False) else 0.0
    spans = np.asarray(sec["span"], dtype=float)
    radii = rh + spans * (rt - rh)
    rmean = 0.5 * (rh + rt)
    angles = free_vortex_angles(rpm, radii, rmean, row["axial_velocity"], row["inlet_angle_mean"], row["exit_angle_mean"])
    rotating = rpm > 0
    a1 = angles["beta1" if rotating else "alpha1"]
    a2 = angles["beta2" if rotating else "alpha2"]
    foils = []
    for i in range(len(spans)):
        ax_c = sec["axial_chord_mm"][i] * 1e-3
        chord = ax_c / np.cos(np.radians(sec["stagger_deg"][i]))
        s_c = (2.0 * np.pi * radii[i] / row["nblades"]) / chord      # local pitch-to-chord ratio
        foils.append(_make_foil(a1[i], a2[i], ax_c, sec["stagger_deg"][i], th, s_c))
    return foils, spans, angles, height, rh, rt, tip_gap


def build_axial_row(case: dict) -> AxialRowCase:
    from microbladegen.aero import Airfoil3D, StackType
    row = case["row"]
    foils, spans, angles, height, rh, rt, tip_gap = axial_row_foils(case)
    rpm = row["rpm"] if row.get("rotating", False) else 0.0
    blade = Airfoil3D(profileArray=foils, profile_loc=list(spans), height=height)
    blade.stack(getattr(StackType, case.get("stack", "Centroid")))
    sw, ln = case.get("sweep"), case.get("lean")
    if sw and any(v != 0 for v in sw["y"]):
        blade.add_sweep(sweep_y=sw["y"], sweep_z=sw["z"])
    if ln and any(v != 0 for v in ln["x"]):
        blade.add_lean(leanX=ln["x"], leanZ=ln["z"])
    b = case.get("build", {})
    blade.build(nProfiles=b.get("nProfiles", 20), num_points=b.get("num_points", 80),
                trailing_edge_points=b.get("trailing_edge_points", 20))

    gas = Gas(gamma=row.get("gamma", 1.4), R=row.get("R", 287.05))
    flow = RowFlow(rpm=rpm, hub_radius=rh, nblades=row["nblades"], axial_velocity=row["axial_velocity"],
                   inlet_angle=row["inlet_angle_mean"], exit_angle=row["exit_angle_mean"],
                   Tt_in=row.get("Tt_in", 288.15), Pt_in=row.get("Pt_in", 101325.0),
                   Tt_out=row.get("Tt_out"), gas=gas)
    metrics = axial_row_metrics(blade, flow)
    return AxialRowCase(case, blade, flow, angles, metrics, row["nblades"], rh, rt, tip_gap)


def build_impeller_case(case: dict, *, tip_clearance: Optional[float] = None) -> ImpellerCase:
    sizing = size_impeller(**case["requirements"])
    geo = case.get("geometry", {})
    cl = case.get("clearance", {})
    budget = clearance_budget(sizing.r2, sizing.u2, material=cl.get("material", "titanium_6al4v"),
                              rotor_delta_T=cl.get("rotor_delta_T", 60.0), casing_delta_T=cl.get("casing_delta_T", 20.0))
    gap = budget.recommended if tip_clearance is None else tip_clearance
    imp = build_impeller(sizing, tip_clearance=gap, tip_thickness=geo.get("tip_thickness_mm", 0.7) * 1e-3,
                         root_thickness_factor=geo.get("root_thickness_factor", 1.5))
    effect = clearance_effect(gap, sizing.b2, "centrifugal")
    return ImpellerCase(case, sizing, imp, budget, effect, impeller_metrics(imp))


def build_case(name_or_case):
    """Build a bundled case by name (or from a case dict / JSON path)."""
    case = load_case(name_or_case) if isinstance(name_or_case, str) else name_or_case
    kind = case["type"]
    if kind == "centrifugal_impeller":
        return build_impeller_case(case)
    if kind == "axial_row":
        return build_axial_row(case)
    raise ValueError(f"unknown case type '{kind}'")
