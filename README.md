# AdditionalDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl)

---

## Documentation

The official documentation is available at:

https://santymax98.github.io/AdditionalDistributions.jl/stable/

---

**AdditionalDistributions.jl** extends the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
ecosystem by providing additional discrete, continuous, and multivariate probability
distributions that are not yet available in the base package.

It implements the core `Distributions.jl` interface (`pdf`, `cdf`, `rand`, etc.) for a growing collection of distributions and includes cumulative probability algorithms for multivariate Gaussian and Student’s *t* models.

---

## 🔹 Key Features

- 📈 **Extensive library of distributions** — zero-inflated, generalized, and heavy-tailed families.
- ⚡ **High-performance QMC algorithms** for multivariate CDFs (`MvGaussian`, `MvTStudent`),  
  significantly faster than [MvNormalCDF.jl](https://github.com/JuliaStats/MvNormalCDF.jl)
  with minimal loss in absolute precision (typically within `1e-5`–`1e-6`).
- 🧮 **Core API compatibility** with `Distributions.jl`.
- 🧩 **Research-oriented architecture**, extensible to new distributional forms.
- 🧠 **Reproducibility-focused testing**, with benchmarks aligned with *mvtnorm* (R, Genz & Bretz 2002).

---

## 🚀 Example

```julia
using AdditionalDistributions

# Univariate and multivariate distributions
d1 = Lomax(α=2.0, λ=3.0)
d2 = ZINB(r=4, θ=0.7, p=0.2)

pdf(d1, 1.5), cdf(d2, 3)

# Accurate CDF for a multivariate t distribution
Σ = [1.0 0.5; 0.5 1.0]
d3 = MvTStudent(ν=5, Σ)
cdf(d3, [-1.0, -1.0], [1.0, 1.0])
```

---

## 🧩 Performance Highlights

| Algorithm                | Library                      | Mean Error | Relative Speed     |
| ------------------------ | ---------------------------- | ---------- | ------------------ |
| QMC–Sobol (this package) | `AdditionalDistributions.jl` | `≈ 1e-5`   | **1.5×–3× faster** |
| Adaptive Genz–Bretz      | `MvNormalCDF.jl`             | `≈ 1e-6`   | slower             |
| QRSVN (R `mvtnorm`)      | Reference                    | `≈ 1e-6`   | —                  |

Our implementation sacrifices a marginal amount of absolute precision
for a substantial speedup in moderate to high dimensions (3–25).

---

## 🧭 Roadmap

* [ ] Add `Generalized Hyperbolic` and `Skew-t` families.
* [ ] Implement flexible parameter fitting (`fit_mle`, `fit_map`).
* [ ] Integrate symbolic representations for documentation.
* [ ] Parallelize with threads
* [ ] Maybe GPU-parallelized QMC backend (planned).

---

## 🤝 Contributing

We welcome contributions!
All code follows the design and testing conventions of
[`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl),
ensuring consistency and interoperability.

### Guidelines

* Open an issue to discuss bugs or ideas.
* Use `@testitem`-based testsets (same as `Distributions.jl`).
* Follow `Documenter.jl` docstring style and type annotations.

Pull requests improving:

* distribution coverage,
* algorithmic efficiency,
* or documentation clarity are particularly encouraged.

---

**Author:** Santiago Jiménez
**License:** MIT
**Repository:** [Santymax98/AdditionalDistributions.jl](https://github.com/Santymax98/AdditionalDistributions.jl)
