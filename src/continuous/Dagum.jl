"""
    Dagum(a, b, p)

A *Dagum* distribution is defined by three parameters: `a`, `b`, and `p`. It is commonly used in economics for modeling income distributions. The probability density function (PDF) of the Dagum distribution is given by:

```math
f(x; a, b, p) = \\frac{a \\cdot p}{x} \\left(\\frac{x}{b}\\right)^{a p} \\left(\\left(\\frac{x}{b}\\right)^a + 1\\right)^{-(p+1)}, \\quad x > 0
```

```julia
Dagum()        # equivalent to Dagum(1, 1, 1)
Dagum(a)    # equivalent to Dagum(a, 1, 1)
Dagum(a, b)    # equivalent to Dagum(a, b, 1)

params(d)        # Get the parameters, i.e. (a, b, p)
```

External links:

* [Dagum distribution on Wikipedia](https://en.wikipedia.org/wiki/Dagum_distribution)
"""
struct Dagum{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T
    b::T
    p::T
    Dagum{T}(a, b, p) where {T<:Real} = new{T}(a, b, p)
end

Dagum(a::Integer, b::Integer, p::Integer; check_args::Bool=true) = Dagum(float(a), float(b), float(p); check_args=check_args)

function Dagum(a::Real, b::Real, p::Real; check_args::Bool=true)
    @check_args Dagum (a, a > zero(a)) (b, b > zero(b)) (p, p > zero(p))
    return Dagum{typeof(a)}(a, b, p)
end

Dagum() = Dagum{Float64}(1, 1, 1)
Dagum(a) = Dagum{Float64}(a, 1, 1)
Dagum(a, b) = Dagum{Float64}(a, b, 1)
@distr_support Dagum 0.0 Inf

# parameters

params(d::Dagum) = (d.a, d.b, d.p)
@inline partype(d::Dagum{T}) where {T<:Real} = T

Base.eltype(::Type{Dagum{T}}) where {T} = T

shape(d::Dagum) = d.a, d.p
scale(d::Dagum) = d.b

#Statistic

function moments(d::Dagum, r::Integer)
    a, b, p = d.a, d.b, d.p
    if r == 0
        return 1
    elseif r < a
        return b^r * SpecialFunctions.gamma(1 - r/a) * SpecialFunctions.gamma(p + r/a)/SpecialFunctions.gamma(p)
    else
        return NaN
    end
end

Statistics.mean(d::Dagum) = moments(d, 1)
 
Statistics.var(d::Dagum) = moments(d, 2) - moments(d, 1)^2

Statistics.std(d::Dagum) = sqrt(Statistics.var(d))

Statistics.median(d::Dagum) = d.b * (2^(1/d.p) - 1)^(-1/d.a)

StatsBase.mode(d::Dagum) = d.b * ((d.a * d.p - 1)/(d.a + 1))^(1/d.a)

#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Dagum, x::Real)
    a, b, p = d.a, d.b, d.p
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0

    return (1 + (x / b)^(-a))^(-p)
end

function Distributions.logpdf(d::Dagum, x::Real)
    a, b, p = d.a, d.b, d.p
    _insupport = insupport(d, x)
    if _insupport
        term_1 = log(a) + (-a * p) * log(b) + log(p) + (-1 + a*p) * log(x)
        term_2 = (-1-p) * log(1 + (x/b)^a)
        return term_1 + term_2
    else
        return -Inf
    end
end


function Distributions.quantile(d::Dagum, q::Real)
    a, b, p = d.a, d.b, d.p
    return b * (q^(-1/p) - 1)^(-1/a)
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Dagum)
    a, b, p = d.a, d.b, d.p
    u = rand(rng)
    return b * (u^(-1/p) - 1)^(-1/a)
end