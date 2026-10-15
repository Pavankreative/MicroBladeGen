"""Clean design summary tables for built cases."""
from __future__ import annotations

import pandas as pd

from microbladegen.micro.cases import ImpellerCase, AxialRowCase
from microbladegen.micro.metrics import impeller_flow_table


def _fmt(df: pd.DataFrame) -> str:
    with pd.option_context("display.float_format", lambda x: f"{x:,.3f}", "display.width", 220,
                           "display.max_columns", 40):
        return df.to_string()


def case_summary(case_obj) -> str:
    """Printable summary of a built case (flow table + section tables + clearance + notes)."""
    c = case_obj.case
    out = [f"{c.get('name', 'case')}: {c.get('description', '')}", ""]
    if isinstance(case_obj, ImpellerCase):
        out += ["MEAN-LINE", _fmt(impeller_flow_table(case_obj.sizing).to_frame()), ""]
        out += ["SECTIONS (hub / mid / tip, measured on the built geometry)", _fmt(case_obj.metrics), ""]
        b, e = case_obj.clearance_budget, case_obj.clearance_effect
        out += ["TIP CLEARANCE",
                f"  cold gap (recommended)  {b.recommended*1e3:.3f} mm   (worst-case stack {b.worst_case*1e3:.3f} mm)",
                f"  disk growth             {b.disk_growth*1e3:.3f} mm   thermal mismatch {b.thermal_mismatch*1e3:.3f} mm",
                f"  gap / exit height       {e.ratio_percent:.2f} %   -> efficiency -{e.delta_efficiency_points:.1f} pts, "
                f"stall margin -{e.delta_stall_margin_points:.1f} pts (first-order)  [{e.status}]", ""]
        notes = list(case_obj.sizing.notes)
    elif isinstance(case_obj, AxialRowCase):
        out += ["FLOW ANGLES (deg from axial, free vortex)", _fmt(case_obj.angles), ""]
        out += ["SECTION METRICS", _fmt(case_obj.metrics), ""]
        notes = []
        if case_obj.tip_gap:
            out += [f"TIP GAP  {case_obj.tip_gap*1e3:.3f} mm (blade height reduced accordingly)", ""]
    else:
        raise TypeError("case_summary expects the result of build_case()")
    notes += list(c.get("notes", []))
    if notes:
        out += ["NOTES"] + [f"  - {n}" for n in notes]
    return "\n".join(out)
