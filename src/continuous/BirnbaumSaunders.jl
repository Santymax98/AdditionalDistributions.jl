"""
    BirnbaumSaunders(μ, α, β)

A *BirnbaumSaunders* distribution, also known as the fatigue life distribution, is used to model life data and times to failure. It is defined by three parameters: the location parameter \\( \\mu \\), the shape parameter \\( \\alpha \\), and the scale parameter \\( \\beta \\). The probability density function (PDF) of the Birnbaum-Saunders distribution is given by:

```math
f(x) = \\frac{\\sqrt{\\frac{x - \\mu}{\\beta}} + \\sqrt{\\frac{\\beta}{x - \\mu}}}{2 \\alpha (x - \\mu)} \\phi\\left( \\frac{\\sqrt{\\frac{x - \\mu}{\\beta}} - \\sqrt{\\frac{\\beta}{x - \\mu}}}{\\alpha} \\right)
```

```julia
BirnbaumSaunders()        # equivalent to BirnbaumSaunders(0, 1, 1)
BirnbaumSaunders(σ)       # equivalent to BirnbaumSaunders(0, α, 1)
BirnbaumSaunders(σ,β)       # equivalent to BirnbaumSaunders(0, α, β)

params(d)        # Get the parameters, i.e. (μ , α, β)
```

External links:

* [Birnbaum Saunders distribution on Wikipedia](https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution)
"""
struct BirnbaumSaunders{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    μ::T
    α::T
    β::T
    BirnbaumSaunders{T}(μ, α, β) where {T<:Real} = new{T}(μ, α, β)
end

BirnbaumSaunders(μ::Integer, α::Integer, β::Integer; check_args::Bool=true) = BirnbaumSaunders(float(μ), float(α), float(β); check_args=check_args)

function BirnbaumSaunders(μ::Real, α::Real, β::Real; check_args::Bool=true)
    @check_args BirnbaumSaunders (α, α > zero(α)) (β, β > zero(β))
    μ, α, β = promote(μ, α, β)
    return BirnbaumSaunders{typeof(α)}(μ, α, β)
end

BirnbaumSaunders() = BirnbaumSaunders{Float64}(0, 1, 1)
BirnbaumSaunders(α) = BirnbaumSaunders{Float64}(0, α, 1)
BirnbaumSaunders(α,β) = BirnbaumSaunders{Float64}(0, α, β)

@distr_support BirnbaumSaunders d.μ Inf


# parameters

params(d::BirnbaumSaunders) = (d.μ, d.α, d.β)
@inline partype(d::BirnbaumSaunders{T}) where {T<:Real} = T

location(d::BirnbaumSaunders) = d.μ
scale(d::BirnbaumSaunders) = d.β
shape(d::BirnbaumSaunders) = d.α

Base.eltype(::Type{BirnbaumSaunders{T}}) where {T} = T

#Statistic


Statistics.mean(d::BirnbaumSaunders) = d.μ + d.β * (1 + (d.α^2/2))

Statistics.var(d::BirnbaumSaunders) = (d.α * d.β)^2 * (1 + (5 * d.α^2) / 4)

StatsBase.skewness(d::BirnbaumSaunders) = (4 * d.α * (11 * d.α^2 + 6))/(5 * d.α^2 + 4)^(3/2)

StatsBase.kurtosis(d::BirnbaumSaunders) = 3 + (6 * d.α^2 * (93 * d.α^2 + 40)) / (5 * d.α^2 + 4)^2
#### evaluate functions CDF, PDF, logPDF an CF

function Distributions.cdf(d::BirnbaumSaunders, x::Real)
    μ, α, β = d.μ, d.α, d.β
    _insupport = insupport(d, x)
    if _insupport
        term = (sqrt((x - μ)/β) - sqrt(β/(x - μ)))/α
        return Distributions.cdf(normal_dist, term)
    else
        return 0.0            
    end
end

function Distributions.pdf(d::BirnbaumSaunders, x::Real)
    μ, α, β = d.μ, d.α, d.β
    _insupport = insupport(d, x)
    if _insupport
        term_1 = (sqrt((x - μ)/β) + sqrt(β/(x - μ)))/(2 * α * (x - μ))
        inner_term = (sqrt((x - μ)/β) - sqrt(β/(x - μ)))/α
        return term_1 * Distributions.pdf(normal_dist, inner_term)
    else
        return 0.0
    end
end

function Distributions.logpdf(d::BirnbaumSaunders, x::Real)
    μ, α, β = d.μ, d.α, d.β
    _insupport = insupport(d, x)
    if _insupport
        term_1 = log((sqrt((x - μ) / β) + sqrt(β / (x - μ))) / (2 * α * (x - μ)))
        inner_term = (sqrt((x - μ) / β) - sqrt(β / (x - μ))) / α
        return term_1 + Distributions.logpdf(normal_dist, inner_term)
    else
        return -Inf
    end
end

function Distributions.quantile(d::BirnbaumSaunders, p::Real)
    μ, α, β = d.μ, d.α, d.β
    Zₚ = Distributions.quantile(normal_dist, p)
    inner_term = (α * Zₚ + sqrt(4 + (α * Zₚ)^2))
    return μ + (β/4) * inner_term^2
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::BirnbaumSaunders)
    μ, α, β = d.μ, d.α, d.β
    X = Distributions.rand(rng, Distributions.Normal(0, α / 2))
    T = β*(1 + 2 * X^2 + 2 * X * sqrt(1 + X^2))
    return μ + T
end