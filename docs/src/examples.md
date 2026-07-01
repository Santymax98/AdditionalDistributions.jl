# Examples

This page collects small, reproducible examples showing how
`AdditionalDistributions.jl` can be used together with the Julia statistics
ecosystem.

The examples are available in the `examples/scripts/` directory of the
repository and can be run with:

```bash
julia --project=examples examples/scripts/01_distributions_api.jl
```

To run all examples:

```bash
for f in examples/scripts/*.jl; do
    echo ""
    echo "=============================="
    echo "Running $f"
    echo "=============================="
    julia --project=examples "$f"
done
```

## Basic `Distributions.jl` interface

File:

```text
examples/scripts/01_distributions_api.jl
```

This example checks that several distributions support the standard
`Distributions.jl` interface, including `pdf`, `logpdf`, `cdf`, `quantile`,
`rand`, and moments when available.

It includes examples such as:

```julia
Lomax(2.5, 1.0)
Burr(2.0, 3.0)
Dagum(2.0, 3.0)
Maxwell(1.5)
Gompertz(1.2, 0.8)
PERT(0.0, 0.6, 1.0)
```

## Heavy-tailed risk models

File:

```text
examples/scripts/02_heavy_tail_risk.jl
```

This example uses positive heavy-tailed distributions to simulate loss models
and estimate empirical risk measures such as Value-at-Risk and Expected
Shortfall.

It compares models such as:

```julia
Lomax(2.2, 1_000.0)
Burr(1.6, 2.5, 1_000.0)
Dagum(2.0, 1_000.0, 1.5)
```

## Discrete count models

File:

```text
examples/scripts/03_count_models.jl
```

This example illustrates count distributions such as zero-inflated and
heavy-tailed discrete models.

Examples include:

```julia
ZIP(2.0, 0.35)
ZINB(5, 0.30, 0.45)
BetaNegBinomial(10, 2.0, 5.0)
Delaporte(2.0, 5.0, 0.45)
Yule(2.5)
Zeta(2.2)
```

Some of these distributions may have infinite or undefined theoretical
moments. In those cases, empirical summaries can be unstable for moderate
sample sizes, which is expected for heavy-tailed count models.

## Multivariate Gaussian rectangular probabilities

File:

```text
examples/scripts/04_mvgaussian_cdf.jl
```

This example computes a rectangular probability for a multivariate Gaussian
distribution:

```julia
res = cdf_result(
    d,
    lower,
    upper;
    rng=rng,
    m=100_000,
    abseps=1e-5,
    releps=1e-5,
)
```

The returned object includes the estimated probability, an error estimate, the
number of evaluations, an algorithm label, and an `inform` flag.

## Multivariate Student-t rectangular probabilities

File:

```text
examples/scripts/05_mvstudent_cdf.jl
```

This example computes a rectangular probability for a multivariate Student-t
distribution using randomized quasi-Monte Carlo integration.

The multivariate Student-t case is handled differently from the Gaussian case
because the common radial scale induces dependence even when the scale matrix is
diagonal.

## `HypothesisTests.jl` workflows

File:

```text
examples/scripts/06_hypothesistests.jl
```

This example shows two workflows.

First, it compares two positive heavy-tailed samples using classical
univariate tests:

```julia
ApproximateTwoSampleKSTest(x, y)
MannWhitneyUTest(x, y)
```

Second, it uses simulated multivariate Gaussian samples together with
multivariate tests from `HypothesisTests.jl`:

```julia
UnequalCovHotellingT2Test(X, Y)
EqualCovHotellingT2Test(X, Y)
BartlettTest(X, Y)
```

This demonstrates how `MvGaussian` samples from `AdditionalDistributions.jl`
can be used directly in downstream statistical testing workflows.

## `Copulas.jl` and `SklarDist`

File:

```text
examples/scripts/07_copulas_sklar.jl
```

This example builds a joint distribution using a Gaussian copula and marginal
distributions from `AdditionalDistributions.jl`:

```julia
C = GaussianCopula(Σ)

marginals = (
    Burr(1.6, 2.5, 1_000.0),
    Dagum(2.0, 1_000.0, 1.5),
)

joint = SklarDist(C, marginals)
```

It then simulates from the joint model and estimates marginal summaries,
rank correlation, and the probability of a joint upper-tail event.

## `Turing.jl` Bayesian model

File:

```text
examples/scripts/08_turing_model.jl
```

This example uses `Maxwell` inside a Bayesian model defined with `Turing.jl`:

```julia
@model function maxwell_model(x)
    σ ~ truncated(Normal(1.5, 0.5), 0.1, 5.0)

    for i in eachindex(x)
        x[i] ~ Maxwell(σ)
    end
end
```

A short NUTS chain recovers the true scale parameter from simulated data,
illustrating compatibility with probabilistic programming workflows.