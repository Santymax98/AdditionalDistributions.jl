# Benchmarks

This page documents the numerical validation used for the multivariate rectangular CDF routines in `AdditionalDistributions.jl`.

The goal is not to claim universal superiority over other implementations. Randomized quasi-Monte Carlo performance depends on dimension, dependence structure, integration region, evaluation budget, randomization, package versions, and hardware. The results below document representative development validation for the current numerical defaults.

Generated benchmark CSV files and the larger cross-language experimental harnesses are development artifacts and are intentionally not committed to the package repository.

## Numerical backend

For correlated Gaussian rectangles, the default floating-point path uses:

1. MVSORT variable reordering;
2. the Genz conditional transformation;
3. a component-by-component (CBC) rank-1 lattice;
4. randomized shifts and tent transformation;
5. cached lattice generating vectors;
6. batched integrand evaluation.

The Gaussian QMC dimension is `d - 1`. Diagonal Gaussian covariance matrices are handled by an exact product shortcut.

For multivariate Student-t rectangles, the same conditional Gaussian construction is augmented by one radial chi-square coordinate from the normal-scale-mixture representation. The QMC dimension is therefore `d`. The radial coordinate and the conditional Gaussian coordinates both use the tent transformation.

For common Student-t degrees of freedom, the radial scaling code also has specialized paths for `ν = 1`, `ν = 2`, and `ν = 4`.

The current default numbers of randomized shifts are:

```julia
MvGaussian:  nshifts = 10
MvTStudent:  nshifts = 8
```

For the default floating-point path, the CBC lattice length is chosen as a prime not exceeding the available per-shift budget and the resulting generating vector is cached. Very small budgets and non-floating scalar types retain a Richtmyer fallback.

Because the lattice length is prime-rounded, the actual number of lattice evaluations can be slightly smaller than the nominal integration budget. `CDFResult.neval` retains the requested-budget convention for API compatibility.

## Committed benchmark scripts

Smaller package benchmarks live in the `benchmark/` directory.

### Basic multivariate Gaussian benchmark

```bash
julia --project=benchmark benchmark/run_mvn_basic.jl
```

### Comparison with `MvNormalCDF.jl`

```bash
julia --project=benchmark benchmark/run_mvn_compare_mvnormalcdf.jl
```

Generated result files should be written under `benchmark/results/` and should not be committed.

## Development validation environment

The tables below were produced on the following single-threaded environment:

- Apple M4;
- 16 GiB RAM;
- Julia 1.12.6;
- Python 3.13.3;
- SciPy 1.18.1;
- NumPy 2.5.2;
- R 4.5.3;
- `mvtnorm` 1.3.3.

Five randomized seeds were used for each stochastic method. The nominal budgets were `10_000` and `100_000`. AdditionalDistributions.jl and R `mvtnorm` were run with zero absolute and relative stopping tolerances so that the nominal budget controlled the run. SciPy was evaluated through its public multivariate CDF interfaces with the corresponding `maxpts` budget. Actual internal evaluation counts may therefore differ slightly across implementations.

Runtime values are medians. Julia allocation measurements are reported only for Julia implementations; allocation numbers are not compared across languages.

## Independent reference probabilities

Accuracy was evaluated against structured deterministic references rather than treating any competing QMC implementation as ground truth.

### Gaussian battery

The Gaussian battery contains 48 cases:

- dimensions `d ∈ {3, 10, 20, 50}`;
- equicorrelation and AR(1) dependence;
- `ρ ∈ {0.3, 0.95}`;
- central, lower-tail, and mixed rectangles.

For equicorrelation, the common-factor representation reduces each multivariate Gaussian rectangle to deterministic one-dimensional quadrature.

For AR(1), the Gaussian Markov property reduces the rectangle probability to a sequence of one-dimensional transition integrals evaluated by deterministic Gauss-Legendre quadrature.

All 48 Gaussian reference cases passed the reference-convergence checks.

### Student-t battery

The Student-t battery contains 72 cases:

- dimensions `d ∈ {3, 10}`;
- equicorrelation and AR(1) dependence;
- `ρ ∈ {0.3, 0.95}`;
- `ν ∈ {1, 4, 10}`;
- central, lower-tail, and mixed rectangles.

The Student-t scale-mixture identity is used to integrate over the radial chi-square variable. Conditional on the radial scale, the inner Gaussian rectangle is evaluated by the deterministic equicorrelation or AR(1) reference described above.

All 72 Student-t reference cases passed the reference-convergence checks.

## Gaussian results

The following table summarizes the complete 48-case battery. `med|err|`, RMSE, and `p90` are based on absolute error against the deterministic reference.

| Nominal budget | Method | med\|err\| | RMSE | p90 | Median time | Julia allocations |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 10,000 | AdditionalDistributions.jl | `8.72e-7` | `6.78e-5` | `7.20e-5` | `5.654 ms` | `0.159 MiB` |
| 10,000 | MvNormalCDF.jl | `3.91e-6` | `1.61e-4` | `2.50e-4` | `6.439 ms` | `0.188 MiB` |
| 10,000 | SciPy | `1.10e-6` | `7.65e-5` | `8.28e-5` | `10.393 ms` | — |
| 10,000 | R `mvtnorm` | `2.06e-6` | `1.09e-4` | `1.56e-4` | `27.000 ms` | — |
| 100,000 | AdditionalDistributions.jl | `1.13e-7` | `1.55e-5` | `2.41e-5` | `56.394 ms` | `0.364 MiB` |
| 100,000 | MvNormalCDF.jl | `9.40e-7` | `3.97e-5` | `5.43e-5` | `64.086 ms` | `1.563 MiB` |
| 100,000 | SciPy | `1.75e-7` | `1.67e-5` | `2.03e-5` | `83.352 ms` | — |
| 100,000 | R `mvtnorm` | `3.90e-7` | `8.73e-5` | `7.90e-5` | `85.500 ms` | — |

Case-wise paired comparisons use the median absolute error across the five seeds for each problem.

| Budget | Comparison | AD wins | Other wins | Unresolved | Median paired error ratio AD/other |
| ---: | --- | ---: | ---: | ---: | ---: |
| 10,000 | vs MvNormalCDF.jl | 45 | 2 | 1 | `0.3292` |
| 10,000 | vs SciPy | 22 | 25 | 1 | `1.1000` |
| 10,000 | vs R `mvtnorm` | 30 | 16 | 2 | `0.6502` |
| 100,000 | vs MvNormalCDF.jl | 43 | 4 | 1 | `0.2712` |
| 100,000 | vs SciPy | 30 | 15 | 3 | `0.6675` |
| 100,000 | vs R `mvtnorm` | 43 | 3 | 2 | `0.2606` |

At the larger budget, AdditionalDistributions.jl has lower paired median error than SciPy in 30 cases versus 15 for SciPy, while also having lower median runtime in this environment. SciPy retains a slightly lower absolute-error `p90` at this budget, so the result should be interpreted as a favorable accuracy/runtime trade-off rather than universal dominance.

## Student-t results

The following table summarizes the complete 72-case Student-t battery.

| Nominal budget | Method | med\|err\| | RMSE | p90 | Median relative error | Median time | Julia allocations |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10,000 | AdditionalDistributions.jl | `8.71e-6` | `1.09e-4` | `1.10e-4` | `9.43e-5` | `4.383 ms` | `0.109 MiB` |
| 10,000 | SciPy | `8.99e-6` | `1.12e-4` | `1.57e-4` | `1.10e-4` | `7.266 ms` | — |
| 10,000 | R `mvtnorm` | `1.27e-5` | `1.06e-4` | `1.17e-4` | `1.53e-4` | `22.500 ms` | — |
| 100,000 | AdditionalDistributions.jl | `7.14e-7` | `1.36e-5` | `1.28e-5` | `7.45e-6` | `44.488 ms` | `0.199 MiB` |
| 100,000 | SciPy | `6.74e-7` | `1.71e-5` | `1.42e-5` | `8.61e-6` | `65.324 ms` | — |
| 100,000 | R `mvtnorm` | `2.81e-6` | `4.87e-5` | `5.82e-5` | `2.73e-5` | `94.000 ms` | — |

Case-wise paired comparisons are:

| Budget | Comparison | AD wins | Other wins | Unresolved | Median paired error ratio AD/other |
| ---: | --- | ---: | ---: | ---: | ---: |
| 10,000 | vs SciPy | 38 | 34 | 0 | `0.9330` |
| 10,000 | vs R `mvtnorm` | 45 | 27 | 0 | `0.7065` |
| 100,000 | vs SciPy | 42 | 30 | 0 | `0.8985` |
| 100,000 | vs R `mvtnorm` | 69 | 3 | 0 | `0.2541` |

The specialized radial paths are visible in runtime by degrees of freedom. At a nominal budget of `100_000`, median AdditionalDistributions.jl runtimes across case-level medians were approximately `21.2 ms` for `ν = 1`, `28.6 ms` for `ν = 4`, and `77.9 ms` for the generic `ν = 10` path. The generic chi-square inverse remains a potential future optimization target, but it is not part of the numerical changes documented here.

## Interpreting QMC error estimates

`CDFResult.error` is a randomized-QMC uncertainty estimate, not a deterministic mathematical error bound.

The validation tables above therefore use the actual absolute or relative discrepancy from an independent deterministic reference whenever the structured reference is available. A small reported QMC error should not by itself be interpreted as proof that the returned probability is equally accurate.

`inform = 1` means that the internal estimated error is above the requested tolerance for the available budget. It does not automatically imply that the probability estimate is unusable.

## Reproducibility checklist

When reporting a multivariate CDF benchmark, include:

- package and language versions;
- operating system and CPU;
- number of execution threads;
- distribution and dimension;
- covariance or scale matrix;
- integration bounds;
- degrees of freedom for Student-t;
- nominal integration budget `m`;
- `nshifts`;
- `abseps` and `releps`;
- random seed;
- runtime methodology;
- allocations when comparable;
- `CDFResult.value`;
- `CDFResult.error`;
- `CDFResult.inform`;
- the source and uncertainty of any external reference probability.

A minimal Gaussian call is:

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 10
ρ = 0.5
Σ = fill(ρ, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill(1.0, d)

dist = MvGaussian(zeros(d), Σ)

res = cdf_result(
    dist,
    lower,
    upper;
    m = 100_000,
    nshifts = 10,
    rng = MersenneTwister(1234),
)
```

For Student-t:

```julia
ν = 4.0
dist = MvTStudent(ν, zeros(d), Σ)

res = cdf_result(
    dist,
    lower,
    upper;
    m = 100_000,
    nshifts = 8,
    rng = MersenneTwister(1234),
)
```

## Algorithm references

The conditional Gaussian transformation follows the approach of Alan Genz, *Numerical Computation of Multivariate Normal Probabilities*, Journal of Computational and Graphical Statistics 1(2), 1992, DOI `10.1080/10618600.1992.10477010`.

The fast component-by-component construction follows the rank-1 lattice methodology described by Dirk Nuyens and Ronald Cools, *Fast Component-by-Component Construction, a Reprise for Different Kernels*, in *Monte Carlo and Quasi-Monte Carlo Methods 2004*, 2006, DOI `10.1007/3-540-31186-6_22`.
