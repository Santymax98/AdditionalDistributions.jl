"""
    GaussKuzmin()

A `Gauss-Kuzmin` distribution is a discrete probability distribution that arises as the limiting distribution of the coefficients in the continued fraction expansion of a random variable uniformly distributed on (0, 1). The probability mass function (PMF) of the Gauss-Kuzmin distribution is given by:

```math
P(X = k) = -\\log_2\\left(1 - \\frac{1}{(1 + k)^2}\\right), \\quad k \\in \\{1, 2, 3, \\dots\\}
```
where:

```julia
GaussKuzmin()        # equivalent to GaussKuzmin()

params(d)        # Get the parameters, i.e. none
```

External link:

* [Gauss Kuzmin distributions on Wikipedia](https://en.wikipedia.org/wiki/Gauss%E2%80%93Kuzmin_distribution)
"""
struct GaussKuzmin <: Distributions.DiscreteUnivariateDistribution end

@distr_support GaussKuzmin 1 Inf

params(d::GaussKuzmin) = ()

# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::GaussKuzmin) = Inf
Statistics.var(d::GaussKuzmin) = Inf
Statistics.median(d::GaussKuzmin) = 2
StatsBase.mode(d::GaussKuzmin) = 1

# Evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(::GaussKuzmin, x::Real)
    x < 1 && return 0.0
    k = floor(Int, x)
    return 1 - log2((k + 2) / (k + 1))
end

function Distributions.pdf(::GaussKuzmin, x::Real)
    insupport(GaussKuzmin(), x) || return 0.0
    k = round(Int, x)
    ε = 1 / (k + 1)^2
    return -log1p(-ε) / log(2)                 # -log2(1-ε)
end

Distributions.logpdf(d::GaussKuzmin, x::Real) = insupport(d, x) ? (log(pdf(d, x))) : -Inf

function Distributions.quantile(::GaussKuzmin, p::Real)
    (0.0 ≤ p ≤ 1.0) || throw(ArgumentError("p ∈ [0,1]"))
    r = (2 - 2^(1 - p)) / (2^(1 - p) - 1)
    return max(1, ceil(Int, r))
end

# Sampling
Distributions.rand(rng::Distributions.AbstractRNG, ::GaussKuzmin) = Distributions.quantile(GaussKuzmin(), rand(rng))