# Change log: Modifications to PyTurbo-Aero

**Originator of the Modifications:** Pavan (see `pyproject.toml`).
**Date of alterations:** 2026-10-15 (initial release 0.1.0).
**Statement of derivation:** MicroBladeGen is a Modification, derived directly from Original Software provided by NASA:
PyTurbo-Aero 1.0.2 (<https://github.com/nasa/pyturbo-aero>), distributed under the NASA Open Source Agreement 1.3.
The originator consents to this characterisation of the alterations as a Modification. See `NOTICE.md`.

## 0.1.0 (2026-10-15)

### Renamed and cleaned
- Package `pyturbo` -> `microbladegen`; distribution `pyturbo-aero` -> `microbladegen`; all imports updated.
- Removed upstream packaging and infrastructure that no longer applied: `setup.py` (it named a different license),
  GitLab CI, `build.sh`, Read the Docs config, the Sphinx `docs/` tree, and about 170 MB of generated tutorial output.
- Removed `tests/heat_pipe_test.py`, which imported a module (`pyturbo.cooling`) that is not in the upstream repository.
- `pyproject.toml` rewritten (Python >= 3.10).

### New: `microbladegen/micro/`
- `impeller.py`: mean-line centrifugal impeller sizing (`size_impeller`) and a builder (`build_impeller`) that drives
  `Centrif` with deterministic camber, thickness in metres, and measured-angle refinement.
- `metrics.py`: axial row metrics (solidity, aspect ratio, t/c, blade speed, inlet/exit Mach, tip blade Mach, chord
  Reynolds number, de Haller, turning), free-vortex angle helper, impeller flow table and section metrics.
- `tip_clearance.py`: disk-growth and thermal-growth estimates, clearance budget, first-order efficiency/stall-margin
  effect, and a geometric tip-gap operation for built impellers.
- `cases.py` + `cases/*.json`: four editable design cases (two impellers, a first-stage axial rotor, a turbine nozzle
  guide vane) with realistic micro-turbojet parameters.
- `summary.py`: printable summary tables. `export_cfd.py`: CFD-ready export (sections, STL patches, flow-path curves,
  manifest with periodic sector angle). `__main__.py`: `python -m microbladegen.micro list | summary | export`.
- `gas.py`: ideal gas with Sutherland viscosity (cold air and hot-gas constants).

### New: examples, docs, tests
- `examples/01_micro_impeller.py`, `02_axial_rotor_lean_sweep.py`, `03_wavy_ngv.py`.
- `docs/01_airfoil_2d.md`, `02_stacking_lean_sweep_wavy.md`, `03_centrifugal_impeller.md`, `04_metrics_clearance_export.md`.
- 43 new tests in `tests/micro/`; the 16 upstream tests are retained and pass under the new package name.

### Behaviours of upstream code found while building this (documented; upstream code not changed)
1. `Centrif` ray-intersection camber ignores `TE_Metal_Angle`. The non-ray mode honours both end angles.
2. `Centrif.tip_clearance` reparameterises the span variable but does not move the blade tip off the shroud. The micro
   builder applies the gap geometrically.
3. `Airfoil3D.add_lean` moves sections axially and `add_sweep` moves them tangentially, the reverse of the usual
   turbomachinery naming. Both take fractions of blade height.
4. The mean camber's first and last few stations are distorted by LE/TE rounding, so blade angles must be probed away
   from the edges.
