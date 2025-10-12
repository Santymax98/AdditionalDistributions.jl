"""
    Burr(c, k, λ)

An *Burr* distribution is defined by three parameters: ``c``, ``k``, and ``\\lambda``, where `c` and `k` are shape parameters and ``\\lambda`` is a scale parameter.

The probability density function (PDF) of the Burr distribution is given by:

```math
f(x; c, k, \\lambda) = \\frac{c k}{\\lambda} \\left(\\frac{x}{\\lambda}\\right)^{c-1} \\left(1 + \\left(\\frac{x}{\\lambda}\\right)^c\\right)^{-(k+1)}
```

```julia
Burr()        # equivalent to Burr(1, 1, 1)
Burr(c)       # equivalent to Burr(c, 1, 1)
Burr(c, k)    # equivalent to Burr(c, k, 1)

params(d)        # Get the parameters, i.e. (c, k, λ)
```

External links:

* [Burr distribution on Wikipedia](https://en.wikipedia.org/wiki/Burr_distribution)
"""
struct Burr{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    c::T
    k::T
    λ::T
    Burr{T}(c, k, λ) where {T<:Real} = new{T}(c, k, λ)
end

Burr(c::Integer, k::Integer, λ::Integer; check_args::Bool=true) = Burr(float(c), float(k), float(λ); check_args=check_args)

function Burr(c::Real, k::Real, λ::Real; check_args::Bool=true)
    @check_args Burr (c, c > zero(c)) (k, k > zero(k)) (λ, λ > zero(λ))
    return Burr{typeof(c)}(c, k, λ)
end

Burr() = Burr{Float64}(1, 1, 1)
Burr(c) = Burr{Float64}(c, 1, 1)
Burr(c, k) = Burr{Float64}(c, k, 1)
@distr_support Burr 0.0 Inf

# parameters

params(d::Burr) = (d.c, d.k, d.λ)
@inline partype(d::Burr{T}) where {T<:Real} = T

Base.eltype(::Type{Burr{T}}) where {T} = T

shape(d::Burr) = d.c, d.k
scale(d::Burr) = d.λ

#Statistic
## aux
function moments(d::Burr, r::Integer)
    c, k, λ = d.c, d.k, d.λ
    if r == 0
        return 1
    elseif 0 < r < c * k 
        return λ^r * SpecialFunctions.gamma((c + r)/c) * SpecialFunctions.gamma(k - r/c)/SpecialFunctions.gamma(k)
    else
        return NaN
    end
end

Statistics.mean(d::Burr) = moments(d, 1)

Statistics.var(d::Burr) = moments(d, 2) - moments(d, 1)^2 

Statistics.std(d::Burr) = sqrt(Statistics.var(d))

Statistics.median(d::Burr) = d.λ * (2^(1/d.k) - 1)^(1/d.c)

StatsBase.mode(d::Burr) = d.λ * ((d.c - 1)/(d.c * d.k + 1))^(1/d.c)

StatsBase.skewness(d::Burr) = (2 * moments(d, 1)^3 - 3 * moments(d, 1) * moments(d, 2) + moments(d, 3))/(moments(d, 2) - moments(d, 1)^2)^(3/2)

StatsBase.kurtosis(d::Burr) =  (-3 * moments(d, 1)^4 + 6 * moments(d, 1)^2 * moments(d, 2) - 4 * moments(d, 1) * moments(d, 3) + moments(d, 4))/(moments(d, 2) - moments(d, 1)^2)^2 - 3

#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Burr, x::Real)
    c, k, λ = d.c, d.k, d.λ
    if !insupport(d, x)
        return 0.0
    end
    z = (x / λ)^c
    F = -LogExpFunctions.expm1(-k * log1p(z))
    return ifelse(F < 1.0, F, prevfloat(1.0))
end


function Distributions.logpdf(d::Burr, x::Real)
    c, k, λ = d.c, d.k, d.λ
    _insupport = insupport(d, x)
    if _insupport
        z = x / λ
        return log(c*k/λ) + (c - 1)*log(z) - (k + 1)*log1p(z^c)
    else
        return -Inf
    end
end

function Distributions.quantile(d::Burr, p::Real)
    c, k, λ = d.c, d.k, d.λ
    (0.0 ≤ p ≤ 1.0) || throw(DomainError(p, "p must be in [0, 1]"))
    if p == 0.0
        return 0.0
    elseif p == 1.0
        return Inf
    end
    inner = LogExpFunctions.expm1(-log1p(-p) / k)
    inner = max(inner, 0.0)
    return λ * inner^(1 / c)
end


## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Burr)
    c, k, λ = d.c, d.k, d.λ
    u = rand(rng)
    return λ * (1/(1 - u)^(1/k) - 1)^(1/c)
end