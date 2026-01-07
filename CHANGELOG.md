# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyaragorn/compare/v0.3.0...HEAD

## [v0.3.0] - 2026-01-08
[v0.3.0]: https://github.com/althonos/pyaragorn/compare/v0.2.0...v0.3.0

### Added
- `threshold_scale` configuration option in `RNAFinder` ([#2](https://github.com/althonos/pyaragorn/pull/2), by [@apcamargo](https://github.com/apcamargo)).
- `__repr__` implementation to `TRNAGene` and `TMRNAGene` (by [@apcamargo](https://github.com/apcamargo)).
- Various properties to `RNAFinder` allowing to access the options given in `__init__`.

### Changed
- Distribute wheels in Limited API mode for Python 3.11 and later.


## [v0.2.0] - 2025-06-02
[v0.2.0]: https://github.com/althonos/pyaragorn/compare/v0.1.0...v0.2.0

### Added
- `permuted` property to check whether a `TMRNAGene` is permuted.
- MyPy type hints for the `pyaragorn.lib` module.

### Changed
- Rename `TRNAGene` and `TMRNAGene` to be more consistent.
- Rename `cds` method of `TMRNAGene` to `orf`.


## [v0.1.0] - 2025-05-27
[v0.1.0]: https://github.com/althonos/pyaragorn/compare/239956f...v0.1.0

Initial release.
