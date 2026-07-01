# Developer Notes

## Test tags

The test suite uses `TestItems.jl`. Common tags are:

- `:continuous`
- `:discrete`
- `:multivariate`
- `:reference`
- `:regression`
- `:slow`

Examples:

```julia
@run_package_tests
@run_package_tests (filter = ti -> :multivariate in ti.tags)
@run_package_tests (filter = ti -> :reference in ti.tags)
```

## Multivariate CDF core

The internal entry point is `mvtcdf`. It prepares limits and covariance matrices, applies MVSORT reordering, and dispatches to either the Gaussian or Student's t randomized QMC path.

Important implementation choices:

- Gaussian diagonal rectangles are evaluated exactly by factorization.
- Gaussian correlated rectangles use folded randomized Richtmyer QMC.
- Student's t rectangles add one chi-square scale coordinate.
- The chi-square coordinate is not folded.
- `inform = 1` should be propagated to users rather than suppressed.

## Adding reference tests

Reference tests should include:

- explicit source of the reference value,
- fixed random seed where randomization is involved,
- tolerances justified by the reference method,
- enough metadata to regenerate the value.

Generated benchmark CSV files should not be committed.
