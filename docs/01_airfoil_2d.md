# 2D profile design (`Airfoil2D`)

An `Airfoil2D` is one blade section in an unrolled plane. The camberline is fixed by the metal angles and
stagger; thickness is then added on each side with Bezier control points, the leading edge is blended
for curvature continuity, and the trailing edge is closed with a semicircle.

```python
from microbladegen.aero import Airfoil2D

f = Airfoil2D(alpha1=53, alpha2=38, axial_chord=0.011, stagger=45)   # mid-span compressor rotor section
f.add_le_thickness(0.04)
f.add_ps_thickness(thicknessArray=[0.05, 0.02, -0.01], expansion_ratio=1.2)
f.add_ss_thickness(thicknessArray=[0.12, 0.10, 0.08, 0.07], camberPercent=0.8, expansion_ratio=1.2)
f.match_le_thickness()
f.te_create(radius=0.25e-3, wedge_ss=3.0, wedge_ps=2.4)
f.add_ss_flow_guidance_2(s_c=0.65, n=10)
f.plot2D()
```

## Inputs

| Input | Meaning | Micro-scale guidance |
|---|---|---|
| `alpha1`, `alpha2` | Inlet and outlet **metal** angles, 0-90 deg from the axial direction | For a rotor use relative angles; for a stator absolute angles. `microbladegen.micro.free_vortex_angles` gives hub/mid/tip flow angles from a velocity triangle (set metal = flow, then add incidence and deviation) |
| `axial_chord` | Axial chord (m) | Chord = `axial_chord / cos(stagger)` |
| `stagger` | Angle from the trailing edge to the leading edge (deg) | Sets the throat together with the camber |
| `add_le_thickness(t)` | Bezier leading-edge thickness control | Larger values give a blunter, more incidence-tolerant nose. Small micro blades usually need a blunt nose for robustness and manufacturability |
| `add_ps_thickness(array)` | Pressure-side control-point thicknesses, spaced by `expansion_ratio` | Negative last value pinches the PS toward the trailing edge |
| `add_ss_thickness(array, camberPercent)` | Suction-side control-point thicknesses; `camberPercent` is where the SS straightens | Front-loaded by default (turbine style); pass `aft_loaded=True` for a compressor style |
| `match_le_thickness()` | Matches second-derivative curvature at the LE | Call after the thickness arrays |
| `te_create(radius, wedge_ss, wedge_ps)` | Semicircular TE of `radius` (m) with wedge angles (deg) | Keep `radius` at or above about 0.2-0.3 mm: thinner edges are hard to make and erode quickly |
| `add_ss_flow_guidance_2(s_c, n)` | Straightens the SS toward the throat. **`s_c` is the pitch-to-chord ratio** | Use the real local `s/c = 1/solidity`, not a default. The case builder computes it from blade count and radius |

## What the thickness numbers mean

The thickness arrays are Bezier *control* values, not final thicknesses. The achieved maximum thickness is
not equal to any input value. Measure it and iterate:

```python
from microbladegen.micro import axial_row_metrics   # t_max_mm and t_over_c per section
```

For the bundled cases, arrays of `[0.12, 0.10, 0.08, 0.07]` (SS) gave `t/c` of 5.7 % at the hub to 9.4 % at the
tip of the 64 krpm rotor (0.73 mm and 1.7 mm), and `[0.16, 0.13, 0.11, 0.10]` gave about 14.5 % on the nozzle
guide vane. Adjust the arrays until `t_max_mm` meets your stress and manufacturing limits.

## Checking the section before building 3D

1. `plot2D()`: check the camber direction and that the SS and PS do not cross.
2. `get_points(n)` returns the SS and PS point arrays (`(n, 2)` each) for exporting a section directly.
3. Confirm the pitch-to-chord ratio you passed to `add_ss_flow_guidance_2` matches the row (blade count, radius).

## Sign and direction

`Airfoil2D` builds with flow top-to-bottom by default (`left_to_right=False`). A compressor row reduces the
relative flow angle (`alpha1 > alpha2` for a rotor, `alpha2 = 0` for an axial-exit stator); a turbine nozzle
increases it (`alpha1 = 0`, `alpha2` large). Always look at the plot of a new section type before relying on
its orientation.
