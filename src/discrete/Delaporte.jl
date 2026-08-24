"""
    Delaporte(λ,α,β)

 
A *Delaporte* distribution is a discrete probability distribution that can be viewed as a compound distribution. It combines a Poisson distribution (with mean `λ`) and a Gamma distribution (with shape parameters `α` and `β`). The probability mass function (PMF) of the Delaporte distribution is given by:

```math
P(X = k) = \\sum_{i=0}^{k} \\frac{\\Gamma(\\alpha + i) \\beta^i \\lambda^{k-i} e^{-\\lambda}}{\\Gamma(\\alpha) i! (k-i)!}, \\quad k \\in \\{0, 1, 2, \\dots\\}
```

```julia
Delaporte()        # equivalent to Delaporte(1, 1, 1)
Delaporte(λ)       # equivalent to Delaporte(λ, 1, 1)
Delaporte(λ, α)    # equivalent to Delaporte(r, α, α)

params(d)        # Get the parameters, i.e. (λ, α, β)
```

External link:

* [Delaporte distribution on Wikipedia](https://en.wikipedia.org/wiki/Delaporte_distribution)
"""
struct Delaporte{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    λ::T
    α::T
    β::T
    Delaporte{T}(λ, α, β) where {T<:Real} = new{T}(λ, α, β)
end

Delaporte(λ::Integer, α::Integer, β::Integer; check_args::Bool=true) = Delaporte(float(λ), float(α), float(β); check_args=check_args)

function Delaporte(λ::Real, α::Real, β::Real; check_args::Bool=true)
    @check_args Delaporte (λ, λ > zero(λ)) (α, zero(α) < α) (β, zero(β) < β)
    λ, α, β = promote(λ, α, β)
    return Delaporte{typeof(λ)}(λ, α, β)
end

Delaporte() = Delaporte(1.0, 1.0, 1.0)
Delaporte(λ) = Delaporte(λ, 1.0, 1.0)
Delaporte(λ, α) = Delaporte(λ, α, α)

@distr_support Delaporte 0 Inf


# parameters

params(d::Delaporte) = (d.λ, d.α, d.β)
partype(::Delaporte{T}) where {T} = T

#Statistic

Statistics.mean(d::Delaporte) = d.λ + (d.α * d.β)

Statistics.var(d::Delaporte) = d.λ + (d.α * d.β)*(1 +d.β)

function StatsBase.mode(d::Delaporte)
    λ, α, β = d.λ, d.α, d.β
    z = (α - 1) * β + λ
    
    if z == floor(z)
        return [z, z + 1]
    else
        return floor(z)
    end
end


#### evaluate functions CDF, PDF, logPDF an CF

@inline _isnonneginteger(x::Real) = isfinite(x) && x >= 0 && x == floor(x)

function _logpdf_delaporte(d::Delaporte, x::Integer)
    λ, α, β = d.λ, d.α, d.β
    log_terms = map(0:x) do i
        SpecialFunctions.loggamma(α + i) + i * log(β) + (x - i) * log(λ) -
        SpecialFunctions.loggamma(α) - SpecialFunctions.logfactorial(i) -
        (α + i) * LogExpFunctions.log1p(β) - SpecialFunctions.logfactorial(x - i)
    end
    return -λ + LogExpFunctions.logsumexp(log_terms)
end

function Distributions.cdf(d::Delaporte, x::Real)
    isnan(float(x)) && return NaN
    x < 0 && return 0.0
    x == Inf && return 1.0
    !isfinite(x) && return 0.0

    kmax = floor(Int, x)
    log_terms = [_logpdf_delaporte(d, k) for k in 0:kmax]
    return exp(LogExpFunctions.logsumexp(log_terms))
end

function Distributions.pdf(d::Delaporte, x::Real)
    _isnonneginteger(x) || return 0.0
    return exp(_logpdf_delaporte(d, floor(Int, x)))
end

function Distributions.logpdf(d::Delaporte, x::Real)
    _isnonneginteger(x) || return -Inf
    return _logpdf_delaporte(d, floor(Int, x))
end

function Distributions.quantile(d::Delaporte, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    p == 0 && return 0
    p == 1 && return Inf

    lo = 0
    hi = 10
    while cdf(d, hi) < p
        hi *= 2
    end

    while lo < hi
        mid = div(lo + hi, 2)
        if cdf(d, mid) < p
            lo = mid + 1
        else
            hi = mid
        end
    end

    return lo
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Delaporte)
    λ, α, β = d.λ, d.α, d.β
    
    Z = rand(rng, Distributions.Gamma(α, β))
    Y = rand(rng, Distributions.Poisson(Z))
    X = rand(rng, Distributions.Poisson(λ))
    
    return Y + X
end