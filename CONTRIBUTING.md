# Contributing to AdditionalDistributions.jl

Thanks for considering a contribution.

## Development setup

```bash
git clone https://github.com/Santymax98/AdditionalDistributions.jl.git
cd AdditionalDistributions.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Test organization

Tests use `TestItems.jl` and `TestItemRunner.jl`.

Common tags:

- `:continuous`
- `:discrete`
- `:multivariate`
- `:reference`
- `:regression`
- `:slow`

Run everything:

```julia
@run_package_tests
```

Run a local subset:

```julia
@run_package_tests (filter = ti -> :multivariate in ti.tags)
@run_package_tests (filter = ti -> :reference in ti.tags)
```

## Adding a univariate distribution

A new univariate distribution should implement the relevant `Distributions.jl` methods:

- constructor with argument checks,
- `params`,
- `minimum` / `maximum` through `@distr_support` when appropriate,
- `pdf`, `logpdf`, `cdf`, `quantile`, `rand`,
- moments when available and mathematically well-defined.

Please add both generic tests and at least one explicit reference-value test.

## Multivariate CDF development

The multivariate core is numerical and randomized. Changes should include:

- seeded tests,
- comparison against known reference cases,
- comments explaining changes in convergence, memory, or runtime,
- benchmarks for representative dimensions and covariance structures.

Avoid changing default values such as `m`, `nshifts`, or `batchsize` without benchmark evidence.

## Reporting numerical cases

When reporting a numerical issue, include all of the following:

- distribution type,
- dimension,
- parameters (`μ`, `Σ`, and `ν` when applicable),
- lower and upper limits,
- `m`, `abseps`, `releps`, `nshifts`, `batchsize`,
- RNG type and seed,
- full `CDFResult`,
- external reference, if available.

Minimal reproducible examples are strongly preferred.
