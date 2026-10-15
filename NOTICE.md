# Notice and attribution

**MicroBladeGen is a Modification of PyTurbo-Aero**, Original Software released by the United States Government as
represented by NASA (developed at NASA Glenn Research Center; developer and citation author Paht Juangphanich).

- Upstream project: <https://github.com/nasa/pyturbo-aero> (version 1.0.2 per its `citation.cff`, released 2022-12-01)
- Distribution terms: NASA Open Source Agreement version 1.3 (NOSA 1.3), as stated in the upstream README. The full
  text is in [`LICENSE`](LICENSE). Under NOSA 1.3 clause 3.A the agreement must accompany every copy of this software.
- Modifications, their dates and their originator are described in [`CHANGES.md`](CHANGES.md), as NOSA 3.C requires.

## What comes from upstream

The core geometry engine is upstream code, renamed but otherwise unchanged in behaviour:
`microbladegen/aero/` (`Airfoil2D`, `Airfoil3D`, `AirfoilWavy`, `Centrif`, `Passage2D`, ...),
`microbladegen/helper/`, `microbladegen/ml/`, the upstream tutorials in `tutorials/` and the upstream tests in `tests/`
(import names updated). Parts of those files may still carry upstream notices; do not remove them (NOSA 3.D).

## What is new

`microbladegen/micro/`, `examples/`, `docs/`, and the tests in `tests/micro/` are new work by the contributor named in
`CHANGES.md`, added under the same agreement.

## No endorsement

This project is **not affiliated with, endorsed by, or sponsored by NASA or any prior contributor** (NOSA 3.E).
NASA's name is used only to identify the origin of the Original Software.

## Licensing questions to confirm with upstream

The upstream repository has no `LICENSE` file; its README cites NOSA 1.3, and its legacy packaging metadata
(`setup.py`) named a different license (GNU GPLv3). This project follows the README and ships the NOSA 1.3 text.
If the licensing matters for your use (commercial distribution, publication, employer policy), confirm with the upstream
maintainers or your legal counsel before relying on it. This notice is not legal advice.

## Export

NOSA 3.J: export of goods or technical data from the United States may require an export license. Nothing in this
project grants one. The micro-turbojet design content is general engineering software; check the rules that apply to
you before sharing detailed propulsion designs across borders.
