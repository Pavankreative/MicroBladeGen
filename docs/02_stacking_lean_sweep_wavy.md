# 3D stacking, lean, sweep and wavy blades (`Airfoil3D`, `AirfoilWavy`)

```python
from microbladegen.aero import Airfoil3D
from microbladegen.helper import StackType

blade = Airfoil3D(profileArray=[hub, mid, tip], profile_loc=[0.0, 0.5, 1.0], height=0.0134)
blade.stack(StackType.Trailing_Edge)
blade.add_lean(leanX=[0, 0.05, 0.12], leanZ=[0, 0.5, 1])        # axial displacement
blade.add_sweep(sweep_y=[0, 0, -0.10], sweep_z=[0, 0.5, 1])     # tangential displacement
blade.build(nProfiles=20, num_points=120, trailing_edge_points=20)
```

## Coordinates

`shft_ss` and `shft_ps` have shape `(n_span, n_points, 3)`: **x** axial, **y** tangential arc length,
**z** radial height above the hub (`z = 0` at the hub, `z = height` at the tip). The absolute radius is
`r_hub + z`. `microbladegen.micro.wrap_to_cylinder` converts to true Cartesian coordinates (`theta = y / r`).

## Stacking

| `StackType` | Stacked on | Typical use |
|---|---|---|
| `Leading_Edge` | leading edges | rotors with a fixed leading-edge radius |
| `Centroid` | section centroids | stators and nozzle guide vanes |
| `Trailing_Edge` | trailing edges | rotors where the exit angle must hold across the span |

## Lean and sweep: names are swapped relative to turbomachinery usage

Measured on the 64 krpm rotor case (blade height 13.4 mm, tip-section centroid shift versus an unleaned blade):

| Call | Input | Tip shift | Direction |
|---|---|---|---|
| `add_lean(leanX=...)` | 0.05 | 0.67 mm | **axial** (x) |
| `add_lean(leanX=...)` | 0.20 | 2.68 mm | **axial** (x) |
| `add_sweep(sweep_y=...)` | -0.20 | -2.68 mm | **tangential** (y) |

So in this library:

- `add_lean(leanX)` moves sections **axially**: this is what is normally called **sweep** in turbomachinery.
- `add_sweep(sweep_y)` moves sections **tangentially**: this is what is normally called **lean**.
- Both are **fractions of blade height** (`0.2 x 13.4 mm = 2.68 mm`), applied along the span at the fractions in
  `leanZ` / `sweep_z` (0 hub, 1 tip) and interpolated smoothly.

Use the shift table above to size the displacement for your own blades. Example 2 prints it.

### Micro-scale considerations

- Sweep and lean displacements are small in absolute terms (about 1-3 mm). Keep the section-to-section shift
  gentle: a kink at a 10 mm span is a stress raiser and hard to machine.
- Tangential lean (the library's `sweep_y`) changes the spanwise pressure distribution and the effective blade
  pitch at each radius; judge its benefit with CFD, because the effect on tip leakage depends on the design.
- Re-check solidity and `t_over_c` after leaning: the section chord is unchanged, but the tangential position changes
  the pitch-wise spacing.

## Wavy blades (`AirfoilWavy`)

`AirfoilWavy` builds the blade with `N` sections, then perturbs leading edge, trailing edge and thickness
sinusoidally with `stretch_thickness_chord(...)` or `stretch_thickness_chord_te(...)` (keeps the trailing edge fixed).
Pass one value per section for each of `SSRatio, PSRatio, LERatio, TERatio, LE_wave_angle, TE_wave_angle`.

```python
t = np.linspace(0, 2*np.pi*3, n_prof)            # 3 waves along the span
blade.stretch_thickness_chord_te(0.03*np.sin(t), 0*t, 0.015*np.sin(t), 0*t, 0*t, 0*t)
```

In this library the effect is mostly on the **suction-side shape**; the leading-edge point itself barely moves.
On the bundled nozzle guide vane (16.5 mm chord, about 2.4 mm maximum thickness) a suction-side ratio of 0.03 moved
the surface by at most 0.085 mm, and 0.3 moved it by 0.85 mm: the displacement is linear in the ratio, about 2.8 mm per
unit ratio for that section. Measure the realised displacement for your own blade (example 3 prints it) and scale the
amplitude to what you can manufacture. For small blades the wavelength (span / number of waves) should stay well above
the machining resolution.
