"""CFD-ready geometry export.

Writes a self-describing folder: wrapped section point files, blade surface STL, hub/shroud
meridional curves, a CSV summary and a ``manifest.json`` with units, blade count, periodic
sector angle, rotation and clearance data. Coordinates are Cartesian with the machine axis
along **x**, units **metres**.

Axial rows are built in an unrolled frame (x axial, y tangential arc length, z radial height
above the hub). The export wraps every section onto its cylinder, ``theta = y / r`` with
``r = r_hub + z``, giving true 3D coordinates. This is exact for the section mid-point and
accurate while the section subtends a small angle.
"""
from __future__ import annotations

import json
import os
from typing import Optional

import numpy as np
import pandas as pd


def wrap_to_cylinder(points: np.ndarray, hub_radius: float) -> np.ndarray:
    """Map unrolled (x, y, z) section points onto cylinders: returns Cartesian (x, Y, Z) with Y = r sin(theta), Z = r cos(theta)."""
    p = np.asarray(points, dtype=float)
    r = hub_radius + p[..., 2]
    th = p[..., 1] / r
    return np.stack([p[..., 0], r * np.sin(th), r * np.cos(th)], axis=-1)   # same y = r sin, z = r cos as Centrif


def _grid_to_triangles(grid: np.ndarray, flip: bool = False) -> np.ndarray:
    """Triangulate a structured (n_span, n_chord, 3) surface into (n_tri, 3, 3)."""
    a, b = grid[:-1, :-1], grid[1:, :-1]
    c, d = grid[1:, 1:], grid[:-1, 1:]
    t1 = np.stack([a, b, c], axis=-2)
    t2 = np.stack([a, c, d], axis=-2)
    tri = np.concatenate([t1.reshape(-1, 3, 3), t2.reshape(-1, 3, 3)], axis=0)
    return tri[:, ::-1, :] if flip else tri


def write_stl(path: str, surfaces: list[np.ndarray], flips: Optional[list[bool]] = None) -> int:
    """Write structured surfaces to a binary STL. Returns the triangle count."""
    from stl import mesh
    flips = flips or [False] * len(surfaces)
    tris = np.concatenate([_grid_to_triangles(g, f) for g, f in zip(surfaces, flips)], axis=0)
    m = mesh.Mesh(np.zeros(tris.shape[0], dtype=mesh.Mesh.dtype))
    m.vectors[:] = tris
    m.save(path)
    return int(tris.shape[0])


def _write_sections(folder: str, prefix: str, ss: np.ndarray, ps: np.ndarray) -> list[str]:
    """One text file per span: SS from LE to TE, then PS from TE back to LE (x y z, metres)."""
    files = []
    for j in range(ss.shape[0]):
        loop = np.vstack([ss[j], ps[j][::-1]])
        # drop exact duplicate consecutive points (shared LE/TE)
        keep = np.r_[True, np.linalg.norm(np.diff(loop, axis=0), axis=1) > 1e-12]
        loop = loop[keep]
        name = f"{prefix}_section_{j:03d}.txt"
        np.savetxt(os.path.join(folder, name), loop, fmt="%.8f")
        files.append(name)
    return files


def _write_curve(path: str, curve: np.ndarray, header: str) -> None:
    np.savetxt(path, curve, delimiter=",", header=header, comments="", fmt="%.8f")


def export_impeller_cfd(imp, outdir: str, *, extra: Optional[dict] = None) -> dict:
    """Export a built :class:`~microbladegen.micro.impeller.MicroImpeller`. Returns the manifest."""
    os.makedirs(os.path.join(outdir, "sections"), exist_ok=True)
    s = imp.sizing
    cen = imp.centrif
    manifest = {
        "type": "centrifugal_impeller", "units": "m", "axis": "x", "rpm": s.rpm,
        "n_full_blades": s.n_full, "n_splitter_blades": s.n_splitter,
        "periodic_sector_deg": 360.0 / s.n_full, "tip_clearance_m": imp.tip_clearance,
        "r1s_m": s.r1s, "r1h_m": s.r1h, "r2_m": s.r2, "b2_m": s.b2,
        "mass_flow_kg_s": s.mass_flow, "Tt1_K": s.Tt1, "Pt1_Pa": s.Pt1,
        "design_pressure_ratio_tt": s.pressure_ratio, "tip_mach": s.tip_mach,
        "notes": [
            "Sections are Cartesian (x, y, z) with the axis along x; one text file per span, SS then PS.",
            "Blade surfaces are open shells (SS and PS); close LE/TE/tip in your CAD or mesher.",
            "Mirror the geometry if the opposite rotation direction is required.",
            "hub_curve.csv and shroud_curve.csv are the meridional (x, r) flow-path contours; the shroud is the casing line, so the blade tip sits tip_clearance below it.",
        ],
    }
    files = {}
    surfaces, flips = [], []
    for label, blade in (("main", cen.mainblade), ("splitter", getattr(cen, "splitterblade", None))):
        if blade is None:
            continue
        files[label] = _write_sections(os.path.join(outdir, "sections"), label, blade.ss_cart_pts, blade.ps_cart_pts)
        surfaces += [blade.ss_cart_pts, blade.ps_cart_pts]
        flips += [False, True]
        write_stl(os.path.join(outdir, f"{label}_blade_surface.stl"), [blade.ss_cart_pts, blade.ps_cart_pts], [False, True])
    _write_curve(os.path.join(outdir, "hub_curve.csv"), imp.hub, "x,r")
    _write_curve(os.path.join(outdir, "shroud_curve.csv"), imp.shroud, "x,r")
    manifest["section_files"] = files
    if extra:
        manifest.update(extra)
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


def export_axial_cfd(row_case, outdir: str, *, extra: Optional[dict] = None) -> dict:
    """Export an :class:`~microbladegen.micro.cases.AxialRowCase` (rotor or stator)."""
    os.makedirs(os.path.join(outdir, "sections"), exist_ok=True)
    b = row_case.blade
    rh = row_case.hub_radius
    ss = wrap_to_cylinder(b.shft_ss, rh)
    ps = wrap_to_cylinder(b.shft_ps, rh)
    files = _write_sections(os.path.join(outdir, "sections"), "blade", ss, ps)
    ntri = write_stl(os.path.join(outdir, "blade_surface.stl"), [ss, ps], [False, True])
    fl = row_case.flow
    r_ps = np.hypot(ps[..., 1], ps[..., 2])
    manifest = {
        "type": "axial_row", "units": "m", "axis": "x", "rotating": fl.rpm > 0, "rpm": fl.rpm,
        "n_blades": row_case.nblades, "periodic_sector_deg": 360.0 / row_case.nblades,
        "hub_radius_m": rh, "tip_radius_m": row_case.tip_radius, "casing_radius_m": row_case.tip_radius,
        "blade_tip_radius_m": float(r_ps.max()), "tip_clearance_m": row_case.tip_gap,
        "axial_velocity_m_s": fl.axial_velocity, "Tt_in_K": fl.Tt_in, "Pt_in_Pa": fl.Pt_in,
        "gamma": fl.gas.gamma, "stl_triangles": ntri, "section_files": {"blade": files},
        "notes": [
            "Sections wrapped onto cylinders: theta = y / r with r = r_hub + z (y treated as arc length).",
            "Blade surface STL contains SS and PS patches only; close the tip and hub in CAD or the mesher.",
            "casing_radius_m is the outer flow-path radius; blade_tip_radius_m is smaller by the tip clearance.",
        ],
    }
    if extra:
        manifest.update(extra)
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


def export_case(case_obj, outdir: str) -> dict:
    """Export any built case (impeller or axial row) and write its summary CSV."""
    from microbladegen.micro.cases import ImpellerCase, AxialRowCase
    from microbladegen.micro.metrics import impeller_flow_table
    os.makedirs(outdir, exist_ok=True)
    extra = {"case": case_obj.case.get("name"), "description": case_obj.case.get("description")}
    if isinstance(case_obj, ImpellerCase):
        man = export_impeller_cfd(case_obj.impeller, outdir, extra=extra)
        impeller_flow_table(case_obj.sizing).to_csv(os.path.join(outdir, "summary_flow.csv"), header=True)
        case_obj.metrics.to_csv(os.path.join(outdir, "summary_sections.csv"))
    elif isinstance(case_obj, AxialRowCase):
        man = export_axial_cfd(case_obj, outdir, extra=extra)
        case_obj.metrics.to_csv(os.path.join(outdir, "summary_sections.csv"))
        case_obj.angles.to_csv(os.path.join(outdir, "summary_angles.csv"), index=False)
    else:
        raise TypeError("export_case expects the result of build_case()")
    return man
