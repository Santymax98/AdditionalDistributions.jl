# Changelog

## Unreleased

### Added

- Direct `cdf_result` support for native `Distributions.AbstractMvNormal` and
  `Distributions.AbstractMvTDist` objects.
- One-sided multivariate CDF convenience methods `cdf(d, x)` for
  `MvGaussian` and `MvTStudent`.
- API-level regression tests covering native/wrapper agreement.

### Changed

- `MvGaussian` now participates directly in the
  `Distributions.AbstractMvNormal` interface and delegates only the primitive
  operations required by that abstraction.
- `MvTStudent` is now explicitly typed around
  `Distributions.AbstractMvTDist`, while preserving the distinction between
  Student-t location/scale parameters and statistical mean/covariance.
- Multivariate documentation now describes native `Distributions.jl`
  interoperability and the type-piracy boundary.

## 0.2.0 - 2026-07-01

### Added

- Structured `CDFResult` return type for multivariate rectangular CDFs.
- `cdf_result` API for `MvGaussian` and `MvTStudent`.
- Folded batch randomized QMC implementation for multivariate Gaussian rectangular probabilities.
- Batch randomized QMC implementation for multivariate Student-t rectangular probabilities.
- `nshifts` and `batchsize` controls for multivariate integration.
- Reference tests for selected univariate continuous and discrete distributions.
- Multivariate tests based on Genz-style reference cases and selected R `mvtnorm::pmvt` values.
- Benchmark scripts comparing `MvGaussian` against `MvNormalCDF.jl`.
- Benchmark documentation and reproducibility guidelines.
- Distribution index page in the documentation.

### Changed

- `MvTStudent` now uses a larger default integration budget than `MvGaussian`.
- Multivariate tests use lowercase tags for easier filtering.
- Documentation now describes `inform` codes, reproducibility requirements, benchmark interpretation, and reporting guidelines.
- README now uses the registered package installation command `Pkg.add("AdditionalDistributions")`.

### Fixed

- `Lomax` quantile endpoint behavior at `p = 1`.

### Notes

- `inform = 1` means that the estimated error did not reach the requested tolerance with the current integration budget. It is not an exception or automatic invalidation of the returned value.
