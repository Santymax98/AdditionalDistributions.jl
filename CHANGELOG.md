# Changelog

## Unreleased

### Added

- Structured `CDFResult` return type for multivariate rectangular CDFs.
- `cdf_result` API for `MvGaussian` and `MvTStudent`.
- Folded batch randomized QMC implementation for multivariate Gaussian rectangular probabilities.
- Batch randomized QMC implementation for multivariate Student's t rectangular probabilities.
- `nshifts` and `batchsize` controls for multivariate integration.
- Reference tests for univariate continuous and discrete distributions.
- Multivariate tests based on Genz-style reference cases and selected R `mvtnorm::pmvt` values.
- Benchmark script comparing `MvGaussian` against `MvNormalCDF.jl`.

### Changed

- `MvTStudent` now uses a larger default integration budget than `MvGaussian`.
- Multivariate tests use lowercase tags for easier filtering.
- Documentation now describes `inform` codes and reproducibility requirements.

### Notes

- `inform = 1` means that the estimated error did not reach the requested tolerance with the current integration budget. It is not an exception or automatic invalidation of the returned value.
