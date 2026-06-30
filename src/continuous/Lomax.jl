"""
    Lomax(α, λ)

A `Lomax` distribution, also known as the Pareto Type II distribution, is commonly used in business, economics, and actuarial science. The probability density function (PDF) of the Lomax distribution is given by:

```math
f(x; \\alpha, \\lambda) = \\frac{\\alpha}{\\lambda} \\left(1 + \\frac{x}{\\lambda}\\right)^{-(\\alpha + 1)}, \\quad x \\geq 0
```
where:

- ``\\alpha`` is a shape parameter
- ``\\lambda`` is a scale parameter

```julia
Lomax()        # equivalent to Lomax(1, 1)
Lomax(α)       # equivalent to Lomax(α, 1)

params(d)      # Get the parameters, i.e. (α, λ)
```

External links:

* [Lomax distribution on Wikipedia](https://en.wikipedia.org/wiki/Lomax_distribution)
"""
struct Lomax{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    α::T
    λ::T
    Lomax{T}(α, λ) where {T<:Real} = new{T}(α, λ)
end

Lomax(α::Integer, λ::Integer; check_args::Bool=true) = Lomax(float(α), float(λ); check_args=check_args)

function Lomax(α::Real, λ::Real; check_args::Bool=true)
    @check_args Lomax (α, α > zero(α)) (λ, λ > zero(λ))
    return Lomax{typeof(α)}(α, λ)
end

Lomax() = Lomax{Float64}(1, 1)
Lomax(α) = Lomax{Float64}(α, 1.0)
@distr_support Lomax 0.0 Inf

# parameters

params(d::Lomax) = (d.α, d.λ)
@inline partype(d::Lomax{T}) where {T<:Real} = T

Base.eltype(::Type{Lomax{T}}) where {T} = T

shape(d::Lomax) = d.α
scale(d::Lomax) = d.λ
#Statistic

Statistics.mean(d::Lomax) = d.α > 1 ? d.λ/(d.α - 1.0) : NaN
function Statistics.var(d::Lomax)
    α, λ = d.α, d.λ
    if α > 2
        return (λ^2 * α) / ((α - 1.0)^2 * (α - 2.0))
    
    elseif 1.0 < α <= 2.0 
        return Inf
    else
        return NaN
    end
end

Statistics.median(d::Lomax) = d.λ * (2^(1/d.α) - 1.0)
StatsBase.mode(d::Lomax) = 0.0
StatsBase.skewness(d::Lomax) = d.α > 3 ? (2.0*(1.0 + d.α)/(d.α - 3)) * sqrt((d.α - 2)/d.α) : NaN 
StatsBase.kurtosis(d::Lomax) = d.α > 4 ? (6 *(d.α^3 + d.α^2 - 6*d.α - 2))/(d.α * (d.α - 3) * (d.α - 4)) : NaN 

#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Lomax, x::Real)
    T = promote_type(typeof(float(x)), typeof(float(d.α)), typeof(float(d.λ)))
    (x < 0) && return zero(T)
    α, λ = T(d.α), T(d.λ)
    xx = T(x)
    return -expm1(-α * log1p(xx / λ))
end

function Distributions.pdf(d::Lomax, x::Real)
    insupport(d, x) || return 0.0
    α, λ = d.α, d.λ
    return (α/λ) * (1 + x/λ)^(-(α + 1))
end

function Distributions.logpdf(d::Lomax, x::Real)
    insupport(d, x) || return -Inf
    α, λ = d.α, d.λ
    return log(α/λ) - (α + 1) * log1p(x / λ)
end

function Distributions.quantile(d::Lomax, p::Real)
    (0 ≤ p ≤ 1) || throw(DomainError(p, "p ∈ [0,1]"))
    T = promote_type(typeof(float(p)), typeof(float(d.α)), typeof(float(d.λ)))
    pp = T(p)
    iszero(pp) && return zero(T)
    isone(pp)  && return oftype(T, Inf)

    α, λ = T(d.α), T(d.λ)
    # λ * ((1 - p)^(-1/α) - 1) → stable: λ * expm1( -log1p(-p)/α )
    return λ * expm1(-log1p(-pp) / α)
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Lomax)
    α, λ = d.α, d.λ
    Y = rand(rng, Distributions.Pareto(α, λ))
    return Y - λ
end