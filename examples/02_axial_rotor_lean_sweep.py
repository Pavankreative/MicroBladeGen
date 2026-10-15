"""First-stage micro axial rotor (0.6 kg/s, 64 krpm) with and without tangential lean / axial sweep.

NOTE on the upstream naming: in this library ``add_sweep(sweep_y=...)`` moves sections
*tangentially* (aerodynamic lean) and ``add_lean(leanX=...)`` moves them *axially*
(aerodynamic sweep). Both are fractions of blade height.

Run:  python examples/02_axial_rotor_lean_sweep.py
"""
import copy
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from microbladegen.micro import load_case, build_case, case_summary

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)

base = load_case("axial_rotor_064krpm")
swept = copy.deepcopy(base)
swept["name"] = "axial_rotor_064krpm_swept"
swept["lean"]["x"] = [0.0, 0.05, 0.12]          # axial displacement, fraction of height (0.12 -> 1.6 mm at tip)
swept["sweep"]["y"] = [0.0, 0.0, -0.10]         # tangential displacement, fraction of height

a, b = build_case(base), build_case(swept)
print(case_summary(a))

height = a.blade.zz[-1]
def centroid_shift(blade, j):
    pts = np.vstack([blade.shft_ss[j], blade.shft_ps[j]])
    return pts[:, :2].mean(axis=0) * 1e3
print("\nTip-section centroid shift due to lean/sweep (mm):",
      np.round(centroid_shift(b.blade, -1) - centroid_shift(a.blade, -1), 3),
      f"(blade height {height*1e3:.1f} mm)")

fig, ax = plt.subplots(figsize=(6, 4))
for row, c, lab in ((a, "tab:blue", "baseline"), (b, "tab:red", "axial sweep + tangential lean")):
    for j in (0, -1):
        loop = np.vstack([row.blade.shft_ss[j], row.blade.shft_ps[j][::-1]]) * 1e3
        ax.plot(loop[:, 0], loop[:, 1], color=c, lw=1, label=lab if j == 0 else None)
ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.legend()
ax.set_title("Hub and tip sections")
plt.tight_layout(); plt.savefig("axial_rotor_lean_sweep.png", dpi=150)
