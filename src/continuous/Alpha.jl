"""
    Alpha(α, β)

An *Alpha* distribution is defined by the following probability density function (PDF):

```math
f(x; \\alpha, \\beta) = \\frac{\\beta}{\\sqrt{2\\pi} \\Phi(\\alpha) x^2} \\exp\\left(-\\frac{(\\alpha - \\frac{\\beta}{x})^2}{2}\\right), \\quad x > 0
```
where:

- ``\\alpha`` is a location parameter
- ``\\beta`` is a scale parameter
- ``\\Phi`` is the cumulative distribution function (CDF) of the standard normal distribution.

```julia
Alpha()        # equivalent to Alpha(1, 1)

params(d)        # Get the parameters, i.e. (α, β)
```
"""
struct Alpha{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    α::T
    β::T
    Alpha{T}(α, β) where {T<:Real} = new{T}(α, β)
end

Alpha(α::Integer, β::Integer; check_args::Bool=true) = Alpha(float(α), float(β); check_args=check_args)

function Alpha(α::Real, β::Real; check_args::Bool=true)
    @check_args Alpha (α, α > zero(α)) (β, β > zero(β))
    α, β = promote(α, β)
    return Alpha{typeof(α)}(α, β)
end

Alpha() = Alpha{Float64}(1, 1)
@distr_support Alpha 0.0 Inf

# parameters

params(d::Alpha) = (d.α, d.β)
@inline partype(d::Alpha{T}) where {T<:Real} = T

Base.eltype(::Type{Alpha{T}}) where {T} = T

location(d::Alpha) = d.α
scale(d::Alpha) = d.β

#Statistic

StatsBase.mode(d::Alpha) = d.β * (sqrt(d.α^2 + 8) - d.α)/4

#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Alpha, x::Real)
    α, β = d.α, d.β
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0
    return Distributions.cdf(normal_dist, α - β / x) / Distributions.cdf(normal_dist, α)
end

function Distributions.pdf(d::Alpha, x::Real)
    α, β = d.α, d.β
    _insupport = insupport(d, x)
    if _insupport
        term_1 = β / (sqrt(2 * π) * Distributions.cdf(normal_dist, α) * x^2)
        term_2 = exp(-(α -(β/x))^2 / 2)
        return term_1 * term_2
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Alpha, x::Real)
    α, β = d.α, d.β
    _insupport = insupport(d, x)
    if _insupport
        term_1 = log(β / (sqrt(2 * π) * Distributions.cdf(normal_dist, α) * x^2))
        term_2 = -(α -(β/x))^2 / 2
        return term_1 + term_2
    else
        return -Inf
    end
end

function Distributions.quantile(d::Alpha, p::Real)
    α, β = d.α, d.β
    x = β / (α - Distributions.quantile(normal_dist, p * Distributions.cdf(normal_dist, α)))
    return x
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Alpha)
    α, β = d.α, d.β
    Y = rand(rng, Distributions.truncated(Distributions.Normal(α, 1); lower=0.0))
    return β / Y
end