"""
    BetaNegBinomial(r,α,β)

A *Beta Negative Binomial* is the compound distribution of the [`NegativeBinomial`](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.NegativeBinomial) distribution where the probability of success `p` is distributed according to the `Beta`. It has three parameters: `r`, the number of successes number of successes until the experiment is stopped and two shape parameters ``\\alpha``, ``\\beta``
 
```math
P(X = k) = \\frac{B(r + k, \\alpha + \\beta)}{B(r, \\alpha) k!} \\frac{\\Gamma(k + \\beta)}{\\Gamma(\\beta)}
```

```julia
BetaNegBinomial()        # equivalent to BetaNegBinomial(1, 1, 1)
BetaNegBinomial(r)       # equivalent to BetaNegBinomial(r, 1, 1)
BetaNegBinomial(r, α)    # equivalent to BetaNegBinomial(r, α, α)

params(d)        # Get the parameters, i.e. (r , α, β)
succprob(d)    # Get the number of successes, i.e. r
```

External links

* [Beta Negative Binomial distribution on Wikipedia](https://en.wikipedia.org/wiki/Beta_negative_binomial_distribution)
"""
struct BetaNegBinomial{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    r::Int
    α::T
    β::T
    BetaNegBinomial{T}(r, α, β) where {T<:Real} = new{T}(r, α, β)
end

BetaNegBinomial(r::Integer, α::Integer, β::Integer; check_args::Bool=true) = BetaNegBinomial(r, float(α), float(β); check_args=check_args)

function BetaNegBinomial(r::Integer, α::Real, β::Real; check_args::Bool=true)
    @check_args BetaNegBinomial (r, r > zero(r)) (α, zero(α) < α) (β, zero(β) < β)
    α, β = promote(α, β)
    return BetaNegBinomial{typeof(α)}(r, α, β)
end

BetaNegBinomial() = BetaNegBinomial(1, 1.0, 1.0)
BetaNegBinomial(r::Integer) = BetaNegBinomial(r, 1.0, 1.0)
BetaNegBinomial(r::Integer, α) = BetaNegBinomial(r, α, α)

@distr_support BetaNegBinomial 0 Inf


# parameters
succprob(d::BetaNegBinomial) = d.r

params(d::BetaNegBinomial) = (d.r, d.α, d.β)
partype(::BetaNegBinomial{T}) where {T} = T

#Statistic


Statistics.mean(d::BetaNegBinomial) = d.α > 1 ? (d.r * d.β)/(d.α - 1) : Inf

function Statistics.var(d::BetaNegBinomial)
    r, α, β = d.r, d.α, d.β
    if α > 2
        numerator = r * β * (r + α - 1) * (β + α - 1)
        denominator = (α - 2) * (α - 1)^2
        return numerator / denominator
    else
        return Inf
    end
end

function StatsBase.skewness(d::BetaNegBinomial)
    r, α, β = d.r, d.α, d.β
    if α > 3
        numerator = (2 * r + α - 1) * (2 * β + α - 1)
        inner_term = (r * β * (r + α - 1) * (β + α - 1))/(α - 2)  
        denominator = (α - 3) * sqrt(inner_term)
        return numerator / denominator
    else
        return Inf
    end
end

#### evaluate functions CDF, PDF, logPDF an CF

@inline _isnonneginteger_bnb(x::Real) = isfinite(x) && x >= 0 && x == floor(x)

function _logpdf_bnb(d::BetaNegBinomial, k::Integer)
    r, α, β = d.r, d.α, d.β
    return SpecialFunctions.loggamma(r + k) + SpecialFunctions.logbeta(α + r, β + k) -
           SpecialFunctions.loggamma(k + 1) - SpecialFunctions.loggamma(r) -
           SpecialFunctions.logbeta(α, β)
end

function Distributions.cdf(d::BetaNegBinomial, x::Real)
    x < 0 && return 0.0
    x == Inf && return 1.0
    !isfinite(x) && return 0.0
    kmax = floor(Int, x)
    log_terms = [_logpdf_bnb(d, k) for k in 0:kmax]
    return exp(LogExpFunctions.logsumexp(log_terms))
end

function Distributions.pdf(d::BetaNegBinomial, x::Real)
    _isnonneginteger_bnb(x) || return 0.0
    return exp(_logpdf_bnb(d, floor(Int, x)))
end

function Distributions.logpdf(d::BetaNegBinomial, x::Real)
    _isnonneginteger_bnb(x) || return -Inf
    return _logpdf_bnb(d, floor(Int, x))
end

function Distributions.quantile(d::BetaNegBinomial, p::Real)
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

function Distributions.cf(d::BetaNegBinomial, t)
    r, α, β = d.r, d.α, d.β
    beta_factor = exp(SpecialFunctions.logbeta(α + r, β) - SpecialFunctions.logbeta(α, β))
    return beta_factor * HypergeometricFunctions._₂F₁(r, β, α + β + r, exp(im * t))
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::BetaNegBinomial)
    r, α, β = d.r, d.α, d.β
    p = rand(rng, Distributions.Beta(α,β))
    return rand(rng, Distributions.NegativeBinomial(r, p))
end