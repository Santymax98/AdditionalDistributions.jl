# ProbabilityDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/ProbabilityDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/ProbabilityDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/ProbabilityDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/ProbabilityDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/ProbabilityDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/ProbabilityDistributions.jl)

---

**ProbabilityDistributions.jl** provides a collection of probability distributions and related functions,  
fully compatible with the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) API.  
It focuses on offering a diverse set of distributions that differ from those established in `Distributions.jl`,  
including accurate multivariate CDF implementations for families such as the **Gaussian** and **Student’s t**,  
which are currently scattered, unimplemented, or not unified across existing libraries.

While compatible with `Distributions.jl` through the same public API (`pdf`, `cdf`, `rand`, etc.),  
**ProbabilityDistributions.jl** is developed as an independent and research-oriented project —  
exploring distributional families and computational methods not currently available in standard Julia libraries.

---

## 🔹 Example

```julia
using ProbabilityDistributions

# Classical and specialized distributions
d1 = Lomax(α=2.0, λ=3.0)
d2 = ZINB(r=4, θ=0.7, p=0.2)

# Evaluate probabilities
pdf(d1, 1.5), cdf(d2, 3)

# Accurate CDF for a multivariate t distribution
Σ = [1.0 0.5; 0.5 1.0]
cdf(MvT(ν=5, Σ), [0.2, 0.3])
```

---

## 🎯 Highlights

* **Wide range of probability distributions**, including zero-inflated and generalized families.
* **Accurate CDFs** for Gaussian, Student’s *t*, and related elliptical models.
* **Numerically stable evaluation** for multivariate and tail probabilities.
* **Compatible API** with `Distributions.jl` (`pdf`, `cdf`, `rand`, etc.).
* **Research-oriented design**, suitable for applied and theoretical modeling.

---

## 🔧 Motivation

Many probability distributions and multivariate CDFs used in applied statistics remain unavailable or fragmented across different Julia libraries.
**ProbabilityDistributions.jl** aims to provide a unified, consistent, and numerically robust framework for these models —
combining clarity, accuracy, and full interoperability within Julia’s statistical ecosystem.

---

*Author:* Santiago Jiménez
*License:* MIT
*Repository:* [Santymax98/ProbabilityDistributions.jl](https://github.com/Santymax98/ProbabilityDistributions.jl)