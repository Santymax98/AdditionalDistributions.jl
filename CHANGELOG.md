# Changelog

## Unreleased

## 0.3.0 - 2026-08-25

### Added

- New continuous distributions: `AsymLaplace`, `HalfCauchy`, `HalfNormal`, and `HalfTDist`.
- New discrete distributions: `PoissonInvGaussian` and `Weibull_Type1`.
- Direct `cdf_result` support for native `Distributions.AbstractMvNormal` and `Distributions.AbstractMvTDist` objects.
- One-sided multivariate CDF methods for `MvGaussian` and `MvTStudent`.
- Fast CBC rank-1 lattice construction and caching for multivariate QMC.

### Changed

- `MvGaussian` now participates directly in the `Distributions.AbstractMvNormal` interface.
- `MvTStudent` is typed around `Distributions.AbstractMvTDist`.
- The default multivariate QMC backend now uses randomized CBC rank-1 lattices with tent transformation.
- Default randomized shifts are now `10` for Gaussian and `8` for Student-t.
- Multivariate numerical documentation and benchmark methodology were expanded.

### Performance

- Added specialized Student-t radial transforms for `ν = 1`, `2`, and `4`.
- Reduced typical runtime and memory usage of multivariate Gaussian and Student-t rectangular CDF evaluation.

### Fixed

- Improved CDF behavior at support boundaries for several univariate distributions.
- Preserved generic non-floating `Real` compatibility in Student-t radial scaling.

## 0.2.1 - 2026-07-01

### Added

- Ecosystem examples for `Distributions.jl`, `HypothesisTests.jl`, `Copulas.jl`, `StatsBase.jl`, and `Turing.jl`.
- `rand(rng, d, n)` and `length(d)` support for `MvGaussian` and `MvTStudent`.

### Changed

- Improved documentation navigation and discoverability.
- Expanded compatibility with newer Julia and dependency versions.

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
