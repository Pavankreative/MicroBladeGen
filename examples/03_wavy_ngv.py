"""Wavy suction-side/edge modification of a micro-turbojet nozzle guide vane.

``AirfoilWavy.stretch_thickness_chord_te`` perturbs the built sections sinusoidally along the span.
In this library the effect is mostly on the suction-side shape: for the bundled vane a suction-side
ratio of 0.03 moves the surface by about 0.085 mm and the displacement is linear in the ratio. At this
scale keep the wavelength (span / waves) well above the machining resolution, and check the realised
displacement (printed below) against what you can manufacture.

Run:  python examples/03_wavy_ngv.py
"""
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from microbladegen.aero import AirfoilWavy
from microbladegen.helper import StackType
from microbladegen.micro import load_case, axial_row_foils

warnings.filterwarnings("ignore")
case = load_case("ngv_micro_turbojet")
foils, spans, angles, height, rh, rt, gap = axial_row_foils(case)

n_prof, n_waves = 60, 3
t = np.linspace(0, 2 * np.pi * n_waves, n_prof)
zero = np.zeros(n_prof)


def build(ss_amp, le_amp):
    b = AirfoilWavy(profileArray=foils, profile_loc=list(spans), height=height)
    b.stack(StackType.Centroid)
    b.build(nProfiles=n_prof, num_points=120, trailing_edge_points=20)
    if ss_amp or le_amp:
        b.stretch_thickness_chord_te(ss_amp * np.sin(t), zero, le_amp * np.sin(t), zero, zero, zero)
    return b


base = build(0.0, 0.0)
wavy = build(0.10, 0.05)

# suction-side displacement of each section relative to the un-waved vane
disp = np.array([np.max(np.linalg.norm(wavy.shft_ss[j] - base.shft_ss[j], axis=1)) for j in range(n_prof)]) * 1e3
span = np.linspace(0, 1, n_prof)
print(f"{n_prof} sections, {n_waves} waves over {height*1e3:.1f} mm "
      f"(wavelength {height/n_waves*1e3:.1f} mm); max suction-side displacement {disp.max():.3f} mm")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
j_pk = int(np.argmax(disp))
for blade, c, lab in ((base, "k", "baseline"), (wavy, "tab:red", f"wavy, span {span[j_pk]:.2f}")):
    j = j_pk if blade is wavy else j_pk
    loop = np.vstack([blade.shft_ss[j], blade.shft_ps[j][::-1]]) * 1e3
    ax[0].plot(loop[:, 0], loop[:, 1], color=c, lw=1, label=lab)
ax[0].set_aspect("equal"); ax[0].legend(); ax[0].set_xlabel("x (mm)"); ax[0].set_title("Section at the peak displacement")
ax[1].plot(span, disp)
ax[1].set_xlabel("span fraction"); ax[1].set_ylabel("max suction-side displacement (mm)")
ax[1].set_title(f"|displacement|: {n_waves} sine waves give {2*n_waves} lobes")
plt.tight_layout(); plt.savefig("wavy_ngv.png", dpi=150)
