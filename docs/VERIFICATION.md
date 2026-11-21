# Verification — 21 November 2026

Python 3.12 on Linux, headless plotting (`MPLBACKEND=Agg`).

- Baseline: 61 tests passed.
- Final: `python -m pytest tests -q` — **75 passed**, 19 warnings (existing SciPy bounded-optimisation tolerance warnings and host pytest/AnyIO warning).
- `python examples/01_micro_impeller.py`, `python examples/02_axial_rotor_lean_sweep.py`, `python examples/03_wavy_ngv.py` — all exit code 0, executed in a separate output directory.

Fixed spline construction that lexicographically sorted points; consecutive duplicate removal now preserves traversal and closed endpoints. Resampling supports two/three-point curves and rejects degenerate/nonfinite data and invalid sample counts. Percentage-spline evaluation now inverts integrated spline arc length and preserves requested sample order. The dependent blade-shell routine trims by path parameter instead of coordinate inequalities; its existing test passes.

Fourteen added regression cases cover query order, straight-line arc positions, repeated points, closed contours, short 2D/3D curves (including three-point cubic splines), degenerate/nonfinite coordinates and invalid counts. Tests establish these behaviours on covered inputs, not universal CAD validity. The full existing geometry suite also passes.

`requirements.txt` delegates dependencies to pyproject.toml. No 3D RANS, structural qualification, experimental engine performance or Windows/8 GB benchmark was run. Target thrust class is application context only. Existing upstream attribution remains applicable.

A legacy three-point cubic-spline special case collapsed two segment lengths into one; removed that obsolete handling and added 2D/3D regressions.
