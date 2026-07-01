# Getting Started

## Installation

AdditionalDistributions.jl is available through Julia's package manager:

```julia
using Pkg
Pkg.add("AdditionalDistributions")
```

Then load it with:

```julia
using AdditionalDistributions
```

Most workflows also use `Distributions.jl`:

```julia
using Distributions
```

## Univariate distributions

AdditionalDistributions.jl provides additional continuous and discrete distributions that follow the familiar `Distributions.jl` interface.

For example, a Lomax distribution:

```julia
using AdditionalDistributions

x = Lomax(2.0, 3.0)

pdf(x, 1.5)
cdf(x, 1.5)
quantile(x, 0.9)
```

A zero-inflated Poisson distribution:

```julia
z = ZIP(2.0, 0.3)

pdf(z, 0)
cdf(z, 4)
rand(z, 10)
```

## Multivariate rectangular probabilities

The package also supports rectangular probabilities for selected multivariate distributions.

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 5
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill( 1.0, d)

dist = MvGaussian(zeros(d), Σ)

res = cdf_result(dist, lower, upper;
    m = 100_000,
    rng = MersenneTwister(1234),
)

res.value
res.error
res.inform
```

Use `cdf(dist, lower, upper)` when you only need the probability estimate. Use `cdf_result` when you need numerical diagnostics.

## Where to go next

- [Distribution Index](Distributions.md)
- [Continuous distributions](bestiary/continuous.md)
- [Discrete distributions](bestiary/discrete.md)
- [Multivariate distributions](bestiary/multivariate.md)
- [Accuracy and Reproducibility](accuracy.md)
- [Benchmarks](benchmarks.md)
