# AdditionalDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl)

**AdditionalDistributions.jl** extends the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) ecosystem with additional continuous, discrete, and multivariate probability distributions.

The package follows the standard `Distributions.jl` interface whenever possible and also provides native Julia routines for rectangular probabilities of multivariate Gaussian and Student-t distributions.

## Installation

AdditionalDistributions.jl is registered in the Julia General registry:

```julia
using Pkg
Pkg.add("AdditionalDistributions")
```

Then:

```julia
using AdditionalDistributions
```

## Highlights

- Additional continuous and discrete probability distributions.
- Zero-inflated, heavy-tailed, reliability, and count models.
- `MvGaussian` integrated with the `Distributions.AbstractMvNormal` interface.
- `MvTStudent` built around native `Distributions.jl` Student-t distributions.
- Rectangular Gaussian and Student-t CDF evaluation through randomized QMC.
- Direct `cdf_result` support for native `MvNormal` and `MvTDist` objects.
- Structured numerical diagnostics through `CDFResult`.
- Reproducible randomized integration with seeded RNGs.
- Integration examples with the wider Julia statistics ecosystem.

## Quick start

```julia
using AdditionalDistributions
using Distributions

d = Lomax(2.0, 3.0)

pdf(d, 1.5)
cdf(d, 1.5)
quantile(d, 0.9)
rand(d)
```

Discrete distributions use the same interface:

```julia
d = ZIP(2.0, 0.3)

pdf(d, 0)
cdf(d, 4)
rand(d)
```

See the [distribution index](https://santymax98.github.io/AdditionalDistributions.jl/stable/Distributions/) for the complete public API.

## Multivariate rectangular probabilities

For a multivariate random vector \(X\), the package evaluates probabilities of the form

```math
P(a_i \le X_i \le b_i,\quad i=1,\ldots,d).
```

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 10
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill(1.0, d)

dist = MvGaussian(zeros(d), Σ)
res = cdf_result(dist, lower, upper; m=100_000, rng=MersenneTwister(1234))

res.value      # probability estimate
res.error      # estimated integration error
res.inform     # integration status
res.neval      # requested integration budget
res.algorithm  # algorithm identifier
```

For package-owned wrappers, `cdf(dist, lower, upper)` returns only the probability estimate.

The same numerical backend can also be used directly with native `Distributions.jl` objects:

```julia
using Distributions

dn = MvNormal(zeros(d), Σ)
dt = MvTDist(4.0, zeros(d), Σ)

cdf_result(dn, lower, upper; rng=MersenneTwister(1234))
cdf_result(dt, lower, upper; rng=MersenneTwister(1234))
```

The package deliberately exposes this interoperability through `cdf_result` rather than adding `Distributions.cdf` methods to external distribution types.

## Numerical methods

Multivariate rectangular probabilities use Genz-style variable conditioning combined with randomized rank-1 lattice quasi-Monte Carlo integration.

The floating-point path uses cached component-by-component lattice construction and tent transformation. Gaussian and Student-t implementations use specialized integration dimensions and Student-t radial transforms where applicable.

Accuracy is assessed against deterministic structured reference probabilities when available. Independent implementations such as `MvNormalCDF.jl`, SciPy, and R's `mvtnorm` are used for comparison rather than as ground truth.

See:

- [Accuracy and reproducibility](https://santymax98.github.io/AdditionalDistributions.jl/stable/accuracy/)
- [Benchmarks](https://santymax98.github.io/AdditionalDistributions.jl/stable/benchmarks/)
- [Multivariate distributions](https://santymax98.github.io/AdditionalDistributions.jl/stable/bestiary/multivariate/)

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Contributing

Contributions are welcome. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for development conventions and guidelines.

## Citation

If you use AdditionalDistributions.jl in your research, please cite:

> S. Jiménez (2025). *AdditionalDistributions.jl — Advanced and Extended Probability Distributions in Julia*.

```bibtex
@misc{Jimenez2025AdditionalDistributions,
    author = {Santiago Jimenez},
    title = {AdditionalDistributions.jl --- Advanced and Extended Probability Distributions in Julia},
    year = {2025},
    url = {https://github.com/Santymax98/AdditionalDistributions.jl},
    note = {Julia package}
}
```

## License

MIT License.
