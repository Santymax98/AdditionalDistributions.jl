# Developer Notes

This page summarizes development conventions for contributors.

## Development setup

```bash
git clone https://github.com/Santymax98/AdditionalDistributions.jl.git
cd AdditionalDistributions.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Test tags

The test suite uses `TestItems.jl` and `TestItemRunner.jl`.

Common tags are:

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

The default test entry point should run the full suite.

## Adding a univariate distribution

A new univariate distribution should implement the relevant `Distributions.jl` methods:

- constructor with parameter checks;
- `params`;
- support through `minimum`, `maximum`, or `@distr_support`;
- `pdf` and `logpdf`;
- `cdf`;
- `quantile`, when available;
- `rand`, when practical;
- moments such as `mean` and `var`, when mathematically well-defined.

Please add:

- generic interface tests;
- at least one explicit reference-value test;
- invalid-parameter tests;
- boundary tests for support endpoints when relevant.

## Multivariate CDF core

The internal multivariate CDF entry point is `mvtcdf`. It prepares limits and covariance matrices, applies MVSORT reordering, and dispatches to either the Gaussian or Student-t randomized QMC path.

Important implementation choices:

- Gaussian diagonal rectangles are evaluated exactly by factorization.
- Gaussian correlated rectangles use folded randomized Richtmyer QMC.
- Student-t rectangles use the normal-scale-mixture representation and add one chi-square coordinate.
- The chi-square coordinate is not folded.
- Diagonal Student-t scale matrices must not be treated as products of independent univariate Student-t probabilities.
- `inform = 1` should be propagated to users rather than suppressed.

## Adding reference tests

Reference tests should include:

- explicit source of the reference value;
- fixed random seed where randomization is involved;
- tolerances justified by the reference method;
- enough metadata to regenerate the value.

For Gaussian rectangular probabilities, Genz-style examples and `MvNormalCDF.jl` comparisons are useful.

For Student-t rectangular probabilities, R's `mvtnorm::pmvt` with `GenzBretz` can be used as an external reference.

## Benchmarks

Benchmark scripts live in `benchmark/`.

Generated benchmark CSV files should not be committed.

When changing numerical defaults such as `m`, `nshifts`, `batchsize`, or the QMC transformation, include benchmark evidence for representative cases.

## Documentation

When adding public API, update:

- docstrings;
- the relevant bestiary page;
- the distribution index if a new type is added;
- tests;
- examples when useful.

Keep README examples short and practical. Put detailed numerical discussion in the documentation pages.
