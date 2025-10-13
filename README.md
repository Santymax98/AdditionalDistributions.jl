# AdditionalDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl)

The [*AdditionalDistributions*](https://github.com/Santymax98/AdditionalDistributions.jl) package is a comprehensive extension of [*Distributions.jl*](https://github.com/JuliaStats/Distributions.jl).
It is designed to expand the functionality of the base package by incorporating both discrete and continuous probability distributions that are not included in `Distributions`, either due to their specialized nature or because they are less commonly used.
*AdditionalDistributions* aims to provide a broader range of statistical tools for data analysis, simulation, and probabilistic modeling, catering to both academic and scientific needs where these additional distributions are essential.

---

## Motivation and Relationship with Distributions.jl

Although several additional distributions have been proposed for inclusion in `Distributions.jl` as early as [JuliaStats/Distributions.jl#124](https://github.com/JuliaStats/Distributions.jl/issues/124) (2013), many of them have not yet been integrated into the main package.
The motivation behind *AdditionalDistributions.jl* is therefore to provide a **flexible yet rigorous, API-compatible environment** for developing and testing both classical and new distribution families — **without compromising statistical quality, numerical stability, or consistency** with the `Distributions.jl` interface.

This package also serves as a foundation for future **multivariate** and **computationally optimized distributions**, designed for efficient integration with [`Copulas.jl`](https://github.com/lrnv/Copulas.jl), where the author is an active co-developer.
These implementations aim to offer faster, more stable probability evaluations for dependence modeling and simulation, while remaining interoperable with the broader JuliaStats ecosystem.

---

## Purpose and Scope

The purpose of `AdditionalDistributions.jl` is to serve as a well-maintained repository for distributions that are too exotic or specialized to be included in the main `Distributions.jl` package.
This includes distributions that are frequently required in specific fields of research but are not yet available in Julia, as well as new or experimental distributions that may not be widely used but have significant potential applications.

For instance, you can define the following distributions (among many others):

```julia
julia> using Distributions, AdditionalDistributions
julia> BetaNegBinomial(r, α, β) # Discrete univariate
julia> Lomax(α, λ)              # Continuous univariate
julia> ZINB(r, θ, p)            # Discrete univariate
julia> Gompertz(η, b)           # Continuous univariate
```

These distributions, along with others like the ARGUS or Zero-Inflated Poisson (ZIP), are essential in various fields of research but are not included in the base `Distributions.jl` package.
`AdditionalDistributions.jl` provides these and more, ensuring that users have access to a broad spectrum of statistical tools.

---

## Key Features

* **Extensive Range of Distributions:** Includes a wide variety of distributions, some of which are not commonly found in standard statistical libraries.
* **Seamless Integration:** Fully compatible with `Distributions.jl`, allowing joint use with other packages in the Julia ecosystem.
* **Detailed Documentation:** Comprehensive references and examples for each implemented distribution.
* **Community-Driven Development:** Contributions are welcome — from new distributions to documentation improvements.

---

## Maintaining Relevance and Utility

Given the importance of these additional distributions in various academic and scientific applications, we are committed to the ongoing maintenance and development of `AdditionalDistributions.jl`.
Our goal is to ensure that the package remains a valuable resource for the Julia community.
To support this, we are considering moving the package under the `JuliaStats` organization, which would provide additional support and visibility, ensuring its long-term sustainability.

With *AdditionalDistributions*, you can:

* **Sample from distributions:** Draw random samples from a variety of distributions.
* **Calculate moments and other properties:** Obtain moments (mean, variance, skewness, kurtosis), entropy, and other statistical metrics.
* **Evaluate probability density/mass functions:** Compute `pdf`, `logpdf`, and related functions.
* **Utilize moment-generating, quantile, and characteristic functions:** Access `mgf`, `quantile`, and `cf` for in-depth statistical analysis.

---

## Future Directions

In future releases, we plan to implement **maximum-likelihood estimators**, expand **multivariate coverage**, and introduce **optimized numerical backends** to support efficient use within `Copulas.jl` and other packages in the JuliaStats and SciML ecosystems.
