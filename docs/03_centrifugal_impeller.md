# Centrifugal impeller design for micro-turbojets

The micro layer wraps the upstream `Centrif` generator with a mean-line sizing step and a builder that makes the
generator's camber inputs behave predictably.

```python
from microbladegen.micro import size_impeller, build_impeller, clearance_budget

s = size_impeller(mass_flow=0.25, rpm=110000, pressure_ratio=3.2, backsweep=30, beta1s_flow=58)
gap = clearance_budget(s.r2, s.u2).recommended
imp = build_impeller(s, tip_clearance=gap)
imp.centrif.mainblade.ss_cart_pts      # (n_span, n_chord, 3) Cartesian points
```

## Mean-line sizing

Inputs: mass flow, speed, impeller total-to-total pressure ratio and efficiency, inlet totals, backsweep, shroud flow
angle, inducer hub-to-shroud ratio, exit flow coefficient `cr2/U2`.

1. **Work.** `dh_is = cp Tt1 (PR^((g-1)/g) - 1)`, `dh0 = dh_is / eta_tt`.
2. **Slip and blade count.** Wiesner slip `sigma = 1 - sqrt(cos(beta2b)) / Z^0.7`. Blade count comes from the Glassman
   heuristic `Z = (pi/30)(110 - alpha2) tan(alpha2)` (alpha2 absolute exit angle from radial), rounded to an even total;
   the two are iterated. With `splitters=True` the total is split into full blades and splitters. Override with `n_full`.
3. **Tip speed.** `dh0 = (1 + f_disk) (sigma - (cr2/U2) tan(beta2b)) U2^2`, with 5 % disk-friction work by default.
   `r2 = U2 / omega`; tip Mach `= U2 / a01`.
4. **Exit width.** Exit static state from total conditions and `dh0`; `b2` from continuity with exit blockage.
5. **Inducer.** With the shroud flow angle fixed, `Cx = U1s / tan(beta1s)`; `r1s` from continuity with inlet blockage.
   The shroud relative Mach and the hub flow angle (constant `Cx`) follow.
6. **Checks.** Diffusion ratio `W1s / W2`, specific speed, machine Reynolds number `rho01 U2 r2 / mu01`, and advisory
   notes (tip Mach above about 1.35, shroud relative Mach above about 1.2, diffusion ratio above 1.9, exit height under
   2.5 mm, tip speed above 480 m/s).

`pressure_ratio` is the **impeller** ratio. Diffuser and volute losses reduce the stage ratio.

Bundled results: the 0.25 kg/s / 110 krpm / PR 3.2 case gives `r2 = 38.0 mm`, `U2 = 438 m/s`, tip Mach 1.29, shroud
relative Mach 0.95, diffusion ratio 1.81, 6 + 6 blades, `b2 = 5.5 mm`.

## Meridional contour

`meridional_endwalls` builds elliptical axial-to-radial turns: the shroud runs from `(0, r1s)` to `(L - b2, r2)` and the
hub (back plate) from `(0, r1h)` to `(L, r2)`, both with a radial tangent at the exit. `L = axial_length_ratio * r2`.

## Blade camber: how the builder drives `Centrif`

The camber is a cubic Bezier in the `(m', theta)` plane, `m' = integral(dm / r)`, whose slope is exactly `tan(beta)`.
Findings that shaped the builder:

- **Use `use_ray_camber=False`.** In ray-intersection mode (used in the upstream tutorials) the trailing-edge metal angle
  has no effect on the result. In the non-ray mode the end tangents equal `LE_Metal_Angle` and `TE_Metal_Angle`.
- **Give both angles the same sign.** A positive LE angle with a negative wrap makes an S-shaped blade.
- **Wrap input.** `wrap_angle` is the angle of the Bezier chord in the `(m', theta)` plane, not the physical wrap. The
  builder sets it to the mean of the inlet and exit angles. The physical wrap comes out at about 70 / 55 / 38 deg at
  hub / mid / tip for the bundled case.
- **Measure away from the edges.** The first and last few stations of the mean camber are distorted by the leading-edge
  and trailing-edge rounding, so angles are probed at stations 8 and -8 (of 100). Refinement passes then nudge the inputs
  until the measured inlet and exit angles match the sizing (within about 1 deg for the bundled case).
- **Thickness units.** `Centrif` takes thickness as a fraction of camber length per side. The builder converts a target tip
  thickness in metres (default 0.7 mm; root 1.5x thicker).

## Tip clearance

Upstream `Centrif.tip_clearance` only reparameterises the span variable: the tip profile is still evaluated on the
shroud, so the blade coordinates do not move. `build_impeller` therefore applies the gap geometrically after the build with
`apply_tip_gap_geometry` (each spanwise column contracted toward the hub so the tip ends `gap` below the casing along
the hub-to-shroud line). The meridional `(mp, theta)` arrays are not updated.

## Limitations

- Mean-line only. There is no loss model for the impeller, no diffuser, no volute and no stress or modal analysis.
- The blade-angle distribution is a single cubic per span; loading distributions that need an S-shaped camber or inducer
  curvature changes need `wrap_displacements` on the upstream profile.
- Splitters start at 50 % of the meridional chord (`splitter_start = 0.5`).
