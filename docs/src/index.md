````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: AdditionalDistributions.jl
  text:
  tagline: An extension of Distributions.jl package. 
  image:
    src: logo.png
    alt: AdditionalDistributions.jl
  actions:
    - theme: brand
      text: Getting started
      link: getting_started
    - theme: alt
      text: View on Github
      link: https://github.com/Santymax98/AdditionalDistributions.jl
    - theme: alt
      text: Bestiary
      link: /bestiary/continuous
---
````
# Welcome to AdditionalDistributions.jl

The [*AdditionalDistributions.jl*](https://github.com/Santymax98/AdditionalDistributions.jl) package is a comprehensive extension of [*Distributions.jl*](https://github.com/JuliaStats/Distributions.jl), designed to expand its functionality by incorporating both **discrete**, **continuous**, and **multivariate** probability distributions not available in the base package.
It provides a unified, research-grade framework for probabilistic modeling, simulation, and numerical evaluation of advanced statistical distributions.

---

## Overview

With *AdditionalDistributions*, you can:

* **Sample from distributions:** Draw random samples from classical and novel models.
* **Compute densities and probabilities:** Access `pdf`, `cdf`, and `logpdf` for a wide range of families.
* **Evaluate multivariate probabilities:** Compute rectangular probabilities for Gaussian and Student’s *t* models using a high-performance QMC integrator.
* **Calculate statistical properties:** Obtain mean, variance, skewness, kurtosis, entropy, and characteristic functions.
* **Extend or define new distributions:** Create custom statistical models compatible with `Distributions.jl`’s interface.

---

## Highlights

| Category                 | Feature                                                | Description                                                    |
| ------------------------ | ------------------------------------------------------ | -------------------------------------------------------------- |
| 🧮 **Compatibility**     | Full API parity with `Distributions.jl`                | Supports `pdf`, `cdf`, `rand`, `entropy`, `params`, etc.       |
| 📈 **Coverage**          | Zero-inflated, generalized, heavy-tailed families      | e.g., `ZINB`, `ZIP`, `Lomax`, `BurrXII`, `BetaNegBinomial`, …  |
| 🔬 **Multivariate CDFs** | QMC-based `cdf` for `MvGaussian` and `MvTStudent`      | Reproducible, accurate results typically within `1e-5`–`1e-6`. |
| 🧩 **Design**            | Modular and research-oriented                          | Mirrors structure, style, and tests from `Distributions.jl`.   |
| 🧠 **Validation**        | Cross-checked against `mvtnorm` (R, Genz & Bretz 2002) | Ensures correctness and reproducibility across ecosystems.     |

---

## Motivation

Many probability distributions and accurate multivariate CDFs remain unavailable or fragmented across Julia packages.
**AdditionalDistributions.jl** aims to unify them under a consistent API, enabling both theoretical and applied work in statistics, econometrics, and data science.

The package is also a testbed for **algorithmic optimization** — including Quasi-Monte-Carlo integration, low-discrepancy sequences, and vectorized numerical transforms — with the goal of eventually contributing mature implementations back to `Distributions.jl`.

---

## Results (illustrative comparison)

Representative results under identical conditions (same machine, same `Σ`, same limits `[a,b]`, same `m`, fixed seed):

| Case                    | Method                        | Probability `p` | Reported error | Note                             |
| ----------------------- | ----------------------------- | --------------: | -------------: | -------------------------------- |
| **MVN** (d=6, m=50 000) | `MvGaussian` (QMC)            |  `0.2130609785` |      `4.19e-5` | `inform = 1` (tolerance not met) |
| **MVN** (d=6, m=50 000) | `MvNormalCDF.jl` (Genz–Bretz) |  `0.2130493052` |      `1.42e-5` | smaller error in this scenario   |

For **Student’s *t*** (e.g., `d=5`, `ν=10`, `m=20 000`), typical outputs are of the form:

```
p ≈ 0.2468,  error ≈ 8.2e-5,  inform = 0
```

**Interpretation**

* QMC yields **reproducible** results (fixed seed) and competitive accuracy.
* Reported errors are **estimators**, not strict bounds; they vary with `m`, dimension, correlation, and region `[a,b]`.
* In some problems Genz–Bretz reports smaller errors for the same `m`; in others QMC matches or surpasses. The goal here is **transparent, honest comparability**.

---

## Usage snippets

### Multivariate Gaussian (rectangular CDF)

```example index
using AdditionalDistributions

μ = zeros(5)
Σ = [1.0 0.6 0.3 0.1 0.0;
     0.6 1.0 0.5 0.3 0.1;
     0.3 0.5 1.0 0.6 0.3;
     0.1 0.3 0.6 1.0 0.6;
     0.0 0.1 0.3 0.6 1.0]

a = fill(-1.0, 5)
b = fill( 1.0, 5)

d = MvGaussian(μ, Σ)
p = cdf(d, a, b; m=50_000)                 # or full=true to get (value, error, inform)
```

### Multivariate Student’s *t* (rectangular CDF)

```example index
ν = 10
Σ = [1.0 0.4 0.2 0.1 0.05;
     0.4 1.0 0.5 0.2 0.1;
     0.2 0.5 1.0 0.4 0.2;
     0.1 0.2 0.4 1.0 0.3;
     0.05 0.1 0.2 0.3 1.0]

a = fill(-0.5, 5)
b = fill( 0.5, 5)

d_t = MvTStudent(ν, zeros(5), Σ)
val, err, info = cdf(d_t, a, b; m=20_000, full=true)
```

---

## Methodological notes

* **Algorithms**:

  * `MvGaussian` / `MvTStudent` use QMC with Richtmyer lattices + Cranley–Patterson shifts; optional antithetic reflection in Gaussian coordinates; reproducibility via seeded RNG.
* **Parameters that matter**:

  * Budget `m`, **dimension** `d`, **correlation** structure, geometry of `[a,b]`. For Student’s *t*, also `ν` and the standardization to correlation.
* **Reporting**:

  * `full=true` returns `(value, error, inform)`. `inform = 1` indicates tolerance not met under the current `m`; increase `m` or relax tolerances.

---

## Roadmap

* Extend coverage to **Generalized Hyperbolic**, **Skew-t**, **Variance-Gamma**.
* Add parameter fitting (`fit_mle`, `fit_map`) for new families.
* Explore GPU-assisted QMC for high-dimensional CDFs.
* Improve diagnostics (`inform`) and optional convergence logging.
* Publish reproducible comparisons across Julia and R ecosystems.

---

## Contributing

We follow the development guidelines and testing conventions of [*Distributions.jl*](https://github.com/JuliaStats/Distributions.jl):

* Use `@testitem` blocks with explicit tolerances and seeded randomness when applicable.
* Include `doctest` examples in docstrings.
* Keep methods type-stable, memory-aware, and numerically robust.

**Ways to contribute**

1. **Report issues** with clear minimal reproductions.
2. **Propose improvements** (performance, modularity, memory) via PRs.
3. **Discuss algorithms** (QMC variance reduction, pivots, transforms).
4. **Documentation** (examples, edge cases, methodological notes).

---

## Related ecosystem

* [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) — Base statistical library.
* [`MvNormalCDF.jl`](https://github.com/JuliaStats/MvNormalCDF.jl) — Genz–Bretz MVN CDF.
* [`StatsBase.jl`](https://github.com/JuliaStats/StatsBase.jl) — Foundational statistics.

---

## Citation

If you use *AdditionalDistributions.jl* in academic work, please cite:

```
S. Jiménez (2025). AdditionalDistributions.jl — Advanced and Extended Probability Distributions in Julia.
Available at: https://github.com/Santymax98/AdditionalDistributions.jl
```

---

**Author:** Santiago Jiménez
**License:** MIT
**Repository:** [Santymax98/AdditionalDistributions.jl](https://github.com/Santymax98/AdditionalDistributions.jl)
