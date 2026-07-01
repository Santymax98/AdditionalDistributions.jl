````@raw html
---
layout: home

title: AdditionalDistributions.jl Documentation

description: Additional probability distributions for Julia compatible with Distributions.jl.

head:
  - - link
    - rel: canonical
      href: https://santymax98.github.io/AdditionalDistributions.jl/stable/

hero:
  name: AdditionalDistributions.jl
  text: Official documentation
  tagline: Additional continuous, discrete, and multivariate probability distributions for Julia.
  image:
    src: logo.png
    alt: AdditionalDistributions.jl
  actions:
    - theme: brand
      text: Getting started
      link: getting_started
    - theme: alt
      text: View on GitHub
      link: https://github.com/Santymax98/AdditionalDistributions.jl
    - theme: alt
      text: Multivariate CDFs
      link: /bestiary/multivariate
---
````

# AdditionalDistributions.jl

`AdditionalDistributions.jl` extends [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) with additional continuous, discrete, and multivariate probability distributions.

The package aims to keep the familiar `Distributions.jl` interface while adding families and numerical routines that are useful for research and applied work.

## Main features

- Extra continuous and discrete univariate distributions.
- Zero-inflated discrete models such as `ZIP`, `ZIB`, and `ZINB`.
- Heavy-tailed distributions such as `Burr`, `Lomax`, and `Zeta`.
- `MvGaussian` rectangular CDF evaluation.
- `MvTStudent` rectangular CDF evaluation.
- Structured multivariate CDF diagnostics through `CDFResult`.

## Multivariate CDFs

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

a = fill(-1.0, d)
b = fill( 1.0, d)

mvnormal = MvGaussian(zeros(d), Σ)
res = cdf_result(mvnormal, a, b; m=100_000, rng=MersenneTwister(1234))

res.value
res.error
res.inform
```

The randomized QMC engine is reproducible when the RNG is fixed. The reported error is an estimate, and `inform = 1` means the requested tolerance was not reached with the current integration budget.

## Documentation map

- [Getting Started](getting_started.md)
- [Compatibility](Compatibility.md)
- [Accuracy and Reproducibility](accuracy.md)
- [Developer Notes](developer.md)
- [Discrete distributions](bestiary/discrete.md)
- [Continuous distributions](bestiary/continuous.md)
- [Multivariate distributions](bestiary/multivariate.md)
