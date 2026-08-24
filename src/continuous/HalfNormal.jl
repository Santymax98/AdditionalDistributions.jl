"""
    HalfNormal(σ)

A `HalfNormal` distribution is the distribution of the absolute value of a
zero-centered normal random variable. Equivalently, if ``Z \\sim N(0, σ^2)``,
then ``X = |Z|`` follows a `HalfNormal(σ)` distribution.

The probability density function (PDF) of the HalfNormal distribution is given by:

```math
f(x; \\sigma) = \\frac{\\sqrt{2}}{\\sigma\\sqrt{\\pi}} \\exp\\left(-\\frac{x^2}{2\\sigma^2}\\right), \\quad x \\geq 0, \\quad \\sigma > 0.
```

```julia
HalfNormal()         # equivalent to HalfNormal(1)
HalfNormal(σ)

params(d)            # Get the parameters, i.e. (σ,)
```

External links:

* [Half-normal distribution on Wikipedia](https://en.wikipedia.org/wiki/Half-normal_distribution)
"""
struct HalfNormal{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    σ::T
    HalfNormal{T}(σ) where {T<:Real} = new{T}(σ)
end

HalfNormal(σ::Integer; check_args::Bool=true) = HalfNormal(float(σ); check_args=check_args)

function HalfNormal(σ::Real; check_args::Bool=true)
    @check_args HalfNormal (σ, σ > zero(σ))
    return HalfNormal{typeof(σ)}(σ)
end

HalfNormal() = HalfNormal(1.0)

@distr_support HalfNormal 0.0 Inf

params(d::HalfNormal) = (d.σ,)
@inline partype(d::HalfNormal{T}) where {T<:Real} = T
Base.eltype(::Type{HalfNormal{T}}) where {T} = T

scale(d::HalfNormal) = d.σ

Statistics.mean(d::HalfNormal) = d.σ * sqrt(2 / π)
Statistics.var(d::HalfNormal) = d.σ^2 * (1 - 2 / π)
Statistics.std(d::HalfNormal) = d.σ * sqrt(1 - 2 / π)
Statistics.median(d::HalfNormal) = d.σ * sqrt(2) * SpecialFunctions.erfinv(0.5)
StatsBase.mode(d::HalfNormal) = zero(partype(d))
StatsBase.skewness(d::HalfNormal) = sqrt(2) * (4 - π) / (π - 2)^(3 / 2)
StatsBase.kurtosis(d::HalfNormal) = 8 * (π - 3) / (π - 2)^2
StatsBase.entropy(d::HalfNormal) = 0.5 * log(2π * ℯ * d.σ^2) - log(2)

function Distributions.cdf(d::HalfNormal, x::Real)
    σ = d.σ
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0

    return SpecialFunctions.erf(x / (sqrt(2) * σ))
end

function Distributions.pdf(d::HalfNormal, x::Real)
    σ = d.σ
    insupport(d, x) || return 0.0
    isinf(x) && return 0.0

    return sqrt(2 / π) / σ * exp(-x^2 / (2 * σ^2))
end

function Distributions.logpdf(d::HalfNormal, x::Real)
    σ = d.σ
    insupport(d, x) || return -Inf
    isinf(x) && return -Inf

    return 0.5 * log(2 / π) - log(σ) - x^2 / (2 * σ^2)
end

function Distributions.quantile(d::HalfNormal, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    p == 0 && return zero(partype(d))
    p == 1 && return oftype(float(d.σ), Inf)

    return d.σ * sqrt(2) * SpecialFunctions.erfinv(p)
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::HalfNormal)
    return abs(rand(rng, Distributions.Normal(0, d.σ)))
end

