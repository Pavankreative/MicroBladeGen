# Micro metrics, tip clearance and CFD export

## Axial row metrics (`axial_row_metrics`)

`RowFlow` describes the mean-line flow through a built `Airfoil3D`; swirl follows a free vortex across the span.

| Column | Definition |
|---|---|
| `solidity` | chord / pitch, pitch = `2 pi r / Z` |
| `aspect_ratio` | blade height / mean chord |
| `t_max_mm`, `t_over_c` | largest nearest-neighbour distance from SS to the PS polyline; over local chord |
| `blade_speed` | `omega r` |
| `inlet_velocity`, `inlet_mach` | `W1` (rotor) or `C1` (stator) over inlet static speed of sound |
| `tip_blade_mach` | `U / a1` (zero for stators); at the tip span it is the blade tip Mach number |
| `reynolds_chord` | `rho W1 c / mu1`, Sutherland viscosity |
| `exit_velocity`, `exit_mach`, `reynolds_chord_exit` | the same at row exit (needs `exit_angle`, `Tt_out`) |
| `velocity_ratio`, `de_haller` | `W2 / W1`; `de_haller` is reported only for diffusing rows (ratio below 1) |
| `turning_deg` | flow turning in the blade frame |

Hot-section rows use `Gas(gamma=1.33)` (see the nozzle guide vane case).

## Impeller metrics

`impeller_flow_table` (mean-line) and `impeller_metrics` (built geometry: blade length, thickness, solidity at the exit
with and without splitters, measured inlet/exit angles, wrap).

## Tip clearance model (`microbladegen.micro.tip_clearance`)

- `disk_rim_growth(radius, tip_speed, material)`: solid-disk rim growth `(1 - nu) rho U^2 r / (4 E)`.
- `thermal_growth(radius, delta_T, material)`: `cte * delta_T * r`.
- `clearance_budget(...)`: stack of manufacturing, runout, out-of-round, disk growth and *differential* thermal growth
  (rotor minus casing). The recommended cold gap is the larger of the RSS stack and the deterministic closure, plus a
  margin.
- `clearance_effect(gap, height, kind)`: first-order efficiency and stall-margin penalties, linear in `gap / height`.
  The defaults are generic literature-range magnitudes (centrifugal: 1.5 efficiency points and 4 stall-margin points per
  1 % of `gap / b2`; axial compressor: 2.0 and 5.0). They give the trend; calibrate them to your hardware.
- `gap_adjusted_height(passage_height, gap)` and `apply_tip_gap_geometry(centrif, gap)` apply a gap to axial and
  centrifugal geometry.

For the bundled 250 g impeller (titanium, `r2 = 38 mm`, `U2 = 438 m/s`) the budget gives a recommended cold gap of about
0.13 mm, of which disk growth is 0.047 mm. That is 2.4 % of the 5.5 mm exit height, which the model flags as `warn`.

## CFD-ready export (`export_case`, `python -m microbladegen.micro export`)

Writes a folder with:

- `sections/*.txt`: one file per span, `x y z` in metres (suction side LE to TE, then pressure side TE to LE), machine axis along **x**.
  Y is `r sin(theta)` and Z is `r cos(theta)`, matching the upstream Cartesian convention.
- `*_blade_surface.stl`: suction- and pressure-side patches (open shells; close the LE, TE and tip in your CAD or mesher).
- `hub_curve.csv`, `shroud_curve.csv` (impellers): meridional `(x, r)` flow-path contours.
- `summary_*.csv`: the flow table and section metrics.
- `manifest.json`: units, rpm, blade counts, **periodic sector angle** (`360 / n_blades`), hub/tip/casing radii, tip
  clearance, and flow conditions.

Axial rows are wrapped onto cylinders with `theta = y / r`, `r = r_hub + z`: exact at the section mid-point and accurate
while a section subtends a small angle. The blade tip sits `tip_clearance` below `casing_radius_m`.

The STL and section files are plain geometry: they do not contain boundary conditions, and the format has not been
tested against a particular mesher.
