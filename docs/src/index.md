````@raw html
---
layout: home

title: AdditionalDistributions.jl Documentation

description: Additional probability distributions for Julia, compatible with Distributions.jl.

head:
  - - link
    - rel: canonical
      href: https://santymax98.github.io/AdditionalDistributions.jl/stable/

hero:
  name: AdditionalDistributions.jl
  text: Additional probability distributions for Julia
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

`AdditionalDistributions.jl` extends the Julia statistics ecosystem with additional probability distributions while following the `Distributions.jl` interface whenever possible.

## Main features

- Continuous and discrete probability distributions beyond the core `Distributions.jl` catalog.
- Zero-inflated, heavy-tailed, reliability, and count models.
- Multivariate Gaussian and Student-t distributions.
- Rectangular multivariate CDF evaluation through randomized QMC.
- Direct `cdf_result` interoperability with native `MvNormal` and `MvTDist`.
- Reproducible numerical integration and structured diagnostics.
- Examples with packages from the Julia statistics ecosystem.

## Installation

```julia
using Pkg
Pkg.add("AdditionalDistributions")
```

## Quick example

```julia
using AdditionalDistributions
using Distributions

d = Lomax(2.0, 3.0)

pdf(d, 1.5)
cdf(d, 1.5)
quantile(d, 0.9)
```

Browse the [Distribution index](Distributions.md) for the available distributions.

## Multivariate probabilities

For a multivariate random vector \(X\), rectangular probabilities have the form

```math
P(a_i \le X_i \le b_i,\quad i=1,\ldots,d).
```

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 5
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill(1.0, d)

dist = MvGaussian(zeros(d), Σ)
res = cdf_result(dist, lower, upper; m=100_000, rng=MersenneTwister(1234))

res.value
res.error
res.inform
```

The same `cdf_result` interface is available for native `Distributions.MvNormal` and `Distributions.MvTDist` objects.

## Numerical validation

The multivariate numerical core is validated with deterministic structured reference probabilities whenever available.

Independent comparisons with `MvNormalCDF.jl`, SciPy, and R's `mvtnorm` are used to evaluate accuracy and performance, but randomized external implementations are not treated as ground truth.

See [Accuracy](accuracy.md), [Benchmarks](benchmarks.md), and [Multivariate distributions](bestiary/multivariate.md) for details.

## Explore

- [Getting Started](getting_started.md)
- [Examples](examples.md)
- [Compatibility](Compatibility.md)
- [Distribution Index](Distributions.md)
- [Accuracy](accuracy.md)
- [Benchmarks](benchmarks.md)
- [Ecosystem integration](ecosystem.md)
- [Developer notes](developer.md)
