# Ecosystem integration

`AdditionalDistributions.jl` is designed to work with the Julia statistics
ecosystem by following the standard `Distributions.jl` interface whenever
possible.

## `Distributions.jl`

The package extends the usual distribution interface:

```julia
pdf(d, x)
logpdf(d, x)
cdf(d, x)
quantile(d, p)
rand(rng, d)
mean(d)
var(d)
```

This makes the distributions usable in many packages that already depend on
the `Distributions.jl` API.

## `HypothesisTests.jl`

Samples generated from distributions in `AdditionalDistributions.jl` can be
used directly with classical statistical tests.

For example, univariate samples can be compared with:

```julia
ApproximateTwoSampleKSTest(x, y)
MannWhitneyUTest(x, y)
```

Multivariate samples from `MvGaussian` can be used with:

```julia
UnequalCovHotellingT2Test(X, Y)
EqualCovHotellingT2Test(X, Y)
BartlettTest(X, Y)
```

See `examples/scripts/06_hypothesistests.jl`.

## `Copulas.jl`

Marginal distributions from `AdditionalDistributions.jl` can be combined with
copulas through `SklarDist`.

```julia
C = GaussianCopula(Σ)

marginals = (
    Burr(1.6, 2.5, 1_000.0),
    Dagum(2.0, 1_000.0, 1.5),
)

joint = SklarDist(C, marginals)
```

See `examples/scripts/07_copulas_sklar.jl`.

## `Turing.jl`

Distributions with compatible `logpdf` methods can be used inside Bayesian
models in `Turing.jl`.

```julia
@model function maxwell_model(x)
    σ ~ truncated(Normal(1.5, 0.5), 0.1, 5.0)

    for i in eachindex(x)
        x[i] ~ Maxwell(σ)
    end
end
```

See `examples/scripts/08_turing_model.jl`.

## Notes

Some advanced distributions may have infinite or undefined moments depending
on their parameter values. This is expected for heavy-tailed models and should
be considered when using empirical summaries, hypothesis tests, or Bayesian
models.