"""
    ZINB(r, θ, p)

Zero-inflated negative binomial distribution. The parameter `r` is the number of successes, `θ` is the structural-zero probability, and `p` is the success probability of the negative binomial component.

```julia
ZINB()       # equivalent to ZINB(1, 0.5, 0.5)
ZINB(r)      # equivalent to ZINB(r, 0.5, 0.5)
ZINB(r, θ)   # equivalent to ZINB(r, θ, 0.5)
params(d)    # returns (r, θ, p)
```
"""
struct ZINB{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    r::Int
    θ::T # structural-zero probability
    p::T # NegativeBinomial success probability
    ZINB{T}(r, θ, p) where {T<:Real} = new{T}(r, θ, p)
end

ZINB(r::Integer, θ::Integer, p::Integer; check_args::Bool=true) =
    ZINB(r, float(θ), float(p); check_args=check_args)

function ZINB(r::Integer, θ::Real, p::Real; check_args::Bool=true)
    θ, p = promote(θ, p)
    @check_args ZINB (r, r > zero(r)) (θ, zero(θ) <= θ <= one(θ)) (p, zero(p) < p <= one(p))
    return ZINB{typeof(θ)}(r, θ, p)
end

ZINB() = ZINB(1, 0.5, 0.5)
ZINB(r) = ZINB(r, 0.5, 0.5)
ZINB(r, θ) = ZINB(r, θ, 0.5)

@distr_support ZINB 0 Inf

params(d::ZINB) = (d.r, d.θ, d.p)
partype(::ZINB{T}) where {T} = T

Statistics.mean(d::ZINB) = (1 - d.θ) * d.r * (1 - d.p) / d.p

function Statistics.var(d::ZINB)
    r, θ, p = d.r, d.θ, d.p
    μ = r * (1 - p) / p
    σ² = r * (1 - p) / p^2
    return (1 - θ) * σ² + θ * (1 - θ) * μ^2
end

@inline _isnonneginteger_zinb(x::Real) = isfinite(x) && x >= 0 && x == floor(x)

function Distributions.cdf(d::ZINB, x::Real)
    isnan(float(x)) && return NaN
    x < 0 && return 0.0
    x == Inf && return 1.0
    !isfinite(x) && return 0.0

    k = floor(Int, x)
    return d.θ + (1 - d.θ) * Distributions.cdf(Distributions.NegativeBinomial(d.r, d.p), k)
end

function Distributions.pdf(d::ZINB, x::Real)
    _isnonneginteger_zinb(x) || return 0.0
    k = floor(Int, x)
    nb = Distributions.NegativeBinomial(d.r, d.p)
    if k == 0
        return d.θ + (1 - d.θ) * Distributions.pdf(nb, 0)
    else
        return (1 - d.θ) * Distributions.pdf(nb, k)
    end
end

function Distributions.logpdf(d::ZINB, x::Real)
    _isnonneginteger_zinb(x) || return -Inf
    k = floor(Int, x)
    nb = Distributions.NegativeBinomial(d.r, d.p)
    if k == 0
        return log(d.θ + (1 - d.θ) * Distributions.pdf(nb, 0))
    else
        return LogExpFunctions.log1p(-d.θ) + Distributions.logpdf(nb, k)
    end
end

function Distributions.quantile(d::ZINB, q::Real)
    (0 <= q <= 1) || throw(DomainError(q, "q must be in [0, 1]"))
    q == 0 && return 0
    q == 1 && return Inf

    lo, hi = 0, 100
    while Distributions.cdf(d, hi) < q
        hi *= 2
    end

    while lo < hi
        mid = (lo + hi) ÷ 2
        if Distributions.cdf(d, mid) < q
            lo = mid + 1
        else
            hi = mid
        end
    end

    return lo
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::ZINB)
    rand(rng) < d.θ && return 0
    return rand(rng, Distributions.NegativeBinomial(d.r, d.p))
end
