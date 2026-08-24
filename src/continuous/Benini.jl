"""
A *Benini* distribution, which is a three-parameter distribution often used in economics and actuarial science. The probability density function (pdf) of the Benini distribution is given by:

```math
f(x; \\alpha, \\beta, \\sigma) = e^{-\\alpha \\log\\left(\\frac{x}{\\sigma}\\right) - \\beta \\left[\\log\\left(\\frac{x}{\\sigma}\\right)\\right]^2} \\cdot \\left(\\frac{\\alpha}{x} + \\frac{2\\beta \\log\\left(\\frac{x}{\\sigma}\\right)}{x}\\right), \\quad x > \\sigma
```

```julia
Benini()        # equivalent to Benini(1, 1, 1)
Benini(α)        # equivalent to Benini(α, 1, 1)
Benini(α, β)        # equivalent to Benini(α, β, 1)

params(d)        # Get the parameters, i.e. (α, β, σ)
```

External links:

* [Benini distribution on Wikipedia](https://en.wikipedia.org/wiki/Benini_distribution)
"""
struct Benini{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    α::T
    β::T
    σ::T
    Benini{T}(α, β, σ) where {T<:Real} = new{T}(α, β, σ)
end

Benini(α::Integer, β::Integer, σ::Integer; check_args::Bool=true) = Benini(float(α), float(β), float(σ); check_args=check_args)

function Benini(α::Real, β::Real, σ::Real; check_args::Bool=true)
    @check_args Benini (α, α > zero(α)) (β, β > zero(β)) (σ, σ > zero(σ))
    return Benini{typeof(α)}(α, β, σ)
end

Benini() = Benini{Float64}(1, 1, 1)
Benini(α) = Benini{Float64}(α, 1, 1)
Benini(α, β) = Benini{Float64}(α, β, 1)
@distr_support Benini d.σ Inf

# parameters

params(d::Benini) = (d.α, d.β, d.σ)
@inline partype(d::Benini{T}) where {T<:Real} = T

Base.eltype(::Type{Benini{T}}) where {T} = T

shape(d::Benini) = (d.α, d.β)
scale(d::Benini) = d.σ
#Statistic


Statistics.mean(d::Benini) = d.σ + (exp(((d.α - 1)^2)/(4*d.β)) * sqrt(π) * d.σ * SpecialFunctions.erfc((d.α - 1)/(2*sqrt(d.β))))/(2*sqrt(d.β))

function Statistics.var(d::Benini)
    α, β, σ = d.α, d.β, d.σ
    term_1 = 4 * exp((α - 2)^2 / (4*β)) * sqrt(β) * SpecialFunctions.erfc((α - 2)/(2 * sqrt(β)))
    term_2 = 4 * exp((α - 1)^2 / (4*β)) * sqrt(β) * SpecialFunctions.erfc((α - 1)/(2 * sqrt(β)))
    term_3 = exp((α - 1)^2 / (2*β)) * sqrt(π) * SpecialFunctions.erfc((α - 1)/(2 * sqrt(β)))^2
    inner_term = term_1 - term_2 - term_3
    return (1/(4 * β)) * sqrt(π) * σ^2 * inner_term
end

Statistics.median(d::Benini) = d.σ * exp((-d.α + sqrt(d.α^2 + d.β * log(16.0)))/(2 * d.β))
#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Benini, x::Real)
    α, β, σ = d.α, d.β, d.σ
    isnan(float(x)) && return NaN
    x <= σ && return 0.0
    isinf(x) && return 1.0

    return 1 - (σ / x)^α * exp(-β * log(x / σ)^2)
end

function Distributions.pdf(d::Benini, x::Real)
    α, β, σ = d.α, d.β, d.σ
    _insupport = insupport(d, x)
    if _insupport
        return (σ/x)^α * exp(-β * (log(x/σ))^2) * (α + 2*β*log(x/σ))/x
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Benini, x::Real)
    α, β, σ = d.α, d.β, d.σ
    _insupport = insupport(d, x)
    if _insupport
        return α * log(σ/x) - β * (log(x/σ))^2 + log(α + 2*β*log(x/σ)) - log(x)
    else
        return -Inf
    end
end

function Distributions.quantile(d::Benini, p::Real)
    α, β, σ = d.α, d.β, d.σ
    term_1 = -α + sqrt(α^2 - 4*β*log(1-p))
    term_2 = 2*β
    return σ * exp(term_1/term_2)
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Benini)
    u = rand(rng)
    return Distributions.quantile(d, u)
end