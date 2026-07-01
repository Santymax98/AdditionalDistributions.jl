````@raw html
---
layout: home

title: AdditionalDistributions.jl Documentation

description: Advanced and extended probability distributions for Julia, compatible with Distributions.jl.

head:
  - - link
    - rel: canonical
      href: https://santymax98.github.io/AdditionalDistributions.jl/stable/

hero:
  name: AdditionalDistributions.jl
  text: Extended probability distributions for Julia
  tagline: Continuous, discrete, and multivariate distributions with a familiar Distributions.jl interface.
  image:
    src: logo.png
    alt: AdditionalDistributions.jl
  actions:
    - theme: brand
      text: Getting started
      link: getting_started
    - theme: alt
      text: Distribution index
      link: Distributions
    - theme: alt
      text: Multivariate CDFs
      link: bestiary/multivariate
    - theme: alt
      text: Benchmarks
      link: benchmarks
---
````

# AdditionalDistributions.jl

`AdditionalDistributions.jl` is a Julia package for additional probability
distributions built on top of the `Distributions.jl` interface.

The package extends the Julia statistics ecosystem with univariate and
multivariate probability distributions that are useful in statistics,
actuarial science, reliability analysis, risk modeling, Bayesian modeling,
simulation, and applied probability.

## Main features

- Additional continuous distributions, including heavy-tailed and reliability models.
- Additional discrete count distributions, including zero-inflated and heavy-tailed models.
- Multivariate Gaussian and Student-t distributions.
- Rectangular CDF utilities for multivariate Gaussian and Student-t distributions.
- Compatibility with the standard `Distributions.jl` API.
- Ecosystem examples with `HypothesisTests.jl`, `Copulas.jl`, and `Turing.jl`.

## Installation

```julia
using Pkg
Pkg.add("AdditionalDistributions")
```

## Documentation overview

Use the guide pages to get started, browse the distribution index, and explore
examples of integration with the Julia statistics ecosystem.

## Quick example

```julia
using AdditionalDistributions
using Distributions

x = Lomax(2.0, 3.0)
pdf(x, 1.5)
cdf(x, 1.5)

z = ZIP(2.0, 0.3)
pdf(z, 0)
cdf(z, 4)
```

See the [Distribution index](Distributions.md) for the public API documented in this package.

## Multivariate rectangular probabilities

For a multivariate distribution, the package evaluates probabilities of rectangles:

```math
P(a_i \leq X_i \leq b_i,\quad i=1,\ldots,d).
```

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 5
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill( 1.0, d)

mvnormal = MvGaussian(zeros(d), Σ)
res = cdf_result(mvnormal, lower, upper;
    m = 100_000,
    rng = MersenneTwister(1234),
)

res.value
res.error
res.inform
```

The randomized QMC engine is reproducible when the RNG is fixed. The reported error is an estimate, and `inform = 1` means the requested tolerance was not reached with the current integration budget.

## Benchmarks and references

Benchmark scripts are available in the repository under `benchmark/`. They compare selected multivariate Gaussian cases against `MvNormalCDF.jl` and provide reproducible performance diagnostics.

Student-t reference values are checked against selected values generated with R's `mvtnorm::pmvt` using the `GenzBretz` algorithm.

See [Benchmarks](benchmarks.md) and [Accuracy and Reproducibility](accuracy.md) for details.

## Citation

If you use AdditionalDistributions.jl in your research, please cite it as:

> S. Jiménez (2025). *AdditionalDistributions.jl — Advanced and Extended Probability Distributions in Julia*. Available at: https://github.com/Santymax98/AdditionalDistributions.jl

BibTeX:

```bibtex
@misc{Jimenez2025AdditionalDistributions,
    author = {Santiago Jimenez},
    title = {AdditionalDistributions.jl --- Advanced and Extended Probability Distributions in Julia},
    year = {2025},
    url = {https://github.com/Santymax98/AdditionalDistributions.jl},
    note = {Julia package}
}
```

If the package has been useful in your work, consider starring the repository. It helps the project gain visibility and supports future development.

## Documentation map

- [Getting Started](getting_started.md)
- [Distribution index](Distributions.md)
- [Compatibility](Compatibility.md)
- [Accuracy and Reproducibility](accuracy.md)
- [Benchmarks](benchmarks.md)
- [Developer Notes](developer.md)
- [Discrete distributions](bestiary/discrete.md)
- [Continuous distributions](bestiary/continuous.md)
- [Multivariate distributions](bestiary/multivariate.md)
