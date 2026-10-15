"""Centrifugal impeller for a ~150-200 N micro-turbojet: size -> build -> check -> export.

Run:  python examples/01_micro_impeller.py [outdir]
"""
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from microbladegen.micro import build_case, case_summary, export_case

warnings.filterwarnings("ignore")
outdir = sys.argv[1] if len(sys.argv) > 1 else "output/impeller_250g"

case = build_case("impeller_250g")          # edit microbladegen/micro/cases/impeller_250g.json to change it
print(case_summary(case))

imp = case.impeller
blade = imp.centrif.mainblade

# Meridional view with the blade tip trace
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].plot(imp.hub[:, 0] * 1e3, imp.hub[:, 1] * 1e3, "k", label="hub")
ax[0].plot(imp.shroud[:, 0] * 1e3, imp.shroud[:, 1] * 1e3, "k--", label="shroud (casing)")
for j, c in ((0, "tab:blue"), (blade.ss_cyl_pts.shape[0] // 2, "tab:green"), (blade.ss_cyl_pts.shape[0] - 1, "tab:red")):
    ax[0].plot(blade.ss_cyl_pts[j, :, 0] * 1e3, blade.ss_cyl_pts[j, :, 1] * 1e3, color=c, lw=1)
ax[0].set_xlabel("x (mm)"); ax[0].set_ylabel("r (mm)"); ax[0].set_aspect("equal"); ax[0].legend()
ax[0].set_title(f"Meridional view, tip gap {imp.tip_clearance*1e3:.2f} mm")

# Blade-angle distribution (mean camber) at hub / mid / tip
from microbladegen.micro.impeller import _mean_camber_angles
for j, n in ((0, "hub"), (blade.ss_cyl_pts.shape[0] // 2, "mid"), (blade.ss_cyl_pts.shape[0] - 1, "tip")):
    ang, _ = _mean_camber_angles(blade, j)
    ax[1].plot(range(len(ang)), -ang, label=n)
ax[1].set_xlabel("chord station"); ax[1].set_ylabel("blade angle magnitude (deg)"); ax[1].legend()
ax[1].set_title("Blade angle: inducer (from axial) to exit (backsweep from radial)")
plt.tight_layout(); plt.savefig("impeller_250g.png", dpi=150)

man = export_case(case, outdir)
print(f"\nCFD package written to {outdir}/ (periodic sector {man['periodic_sector_deg']:.1f} deg)")
