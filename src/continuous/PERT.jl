"""
    PERT(a, b, m)

A `PERT` distribution is commonly used in project management for modeling the uncertainty of task durations. The probability density function (PDF) of the PERT distribution is given by:

```math
f(x; a, b, c) = \\frac{(x - a)^{\\alpha - 1} (c - x)^{\\beta - 1}}{B(\\alpha, \\beta) (c - a)^{\\alpha + \\beta - 1}}, \\quad a \\leq x \\leq c
```
where:

- ``\\alpha = 1 + 4\\frac{b - a}{c - a}``
- ``\\beta = 1 + 4\\frac{c - b}{c - a}``
- ``B(\\alpha, \\beta)`` is a [beta function](https://en.wikipedia.org/wiki/Beta_function)

```julia
PERT()        # equivalent to PERT(0, 0.5, 1)

params(d)        # Get the parameters, i.e. (a, b, m)
```

External links:

* [PERT distribution on Wikipedia](https://en.wikipedia.org/wiki/PERT_distribution)
"""
struct PERT{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T # minimum
    b::T # most likely
    m::T # maximum
    PERT{T}(a, b, m) where {T<:Real} = new{T}(a, b, m)
end

# Constructor for integer inputs
PERT(a::Integer, b::Integer, m::Integer; check_args::Bool=true) = PERT(float(a), float(b), float(m); check_args=check_args)

# Main constructor with argument checking
function PERT(a::Real, b::Real, m::Real; check_args::Bool=true)
    @check_args PERT (b, b > a) (m, m > b)
    return PERT{typeof(a)}(a, b, m)
end

PERT() = PERT{Float64}(0, 0.5, 1)
@distr_support PERT d.a d.m

# Parameters
params(d::PERT) = (d.a, d.b, d.m)
@inline partype(d::PERT{T}) where {T<:Real} = T

Base.eltype(::Type{PERT{T}}) where {T} = T
# Accessors for individual parameters
_α(d::PERT) = 1 + 4 * (d.b - d.a) / (d.m - d.a)
_β(d::PERT) = 1 + 4 * (d.m - d.b) / (d.m - d.a)
# Location and scale (mean and variance)
Statistics.mean(d::PERT) = (d.a + 4d.b + d.m) / 6
Statistics.var(d::PERT) = (Statistics.mean(d) - d.a) * (d.m - Statistics.mean(d)) / 7.0
Statistics.median(d::PERT) = Distributions.quantile(d, 0.5)
# Mode
StatsBase.mode(d::PERT) = d.b
StatsBase.skewness(d::PERT) = (2 * (_β(d) - _α(d)) * sqrt(_α(d) + _β(d) + 1)) / ((_α(d) + _β(d) + 2) * sqrt(_α(d) * _β(d)))
# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::PERT, x::Real)
    a, b, m = d.a, d.b, d.m
    x < a && return 0.0
    x >= m && return 1.0
    z = (x - a) / (m - a)
    alpha, beta = _α(d), _β(d)
    value = SpecialFunctions.beta_inc(alpha, beta, z)
    return value[1]
end

function Distributions.pdf(d::PERT, x::Real)
    a, b, m = d.a, d.b, d.m
    _insupport = insupport(d, x)
    if _insupport
        alpha, beta = _α(d), _β(d)
        num = (x - a)^(alpha - 1) * (m - x)^(beta - 1)
        dem = SpecialFunctions.beta(alpha, beta) * (m - a)^(alpha + beta - 1)
        return num/dem
    else
        return 0.0            
    end
end

function Distributions.logpdf(d::PERT, x::Real)
    a, b, m = d.a, d.b, d.m
    _insupport = insupport(d, x)
    if _insupport
        alpha, beta = _α(d), _β(d)
        term_1 = (alpha - 1) * log(x - a) + (beta - 1) * log(m - x)
        term_2 = SpecialFunctions.logbeta(alpha, beta) + (alpha + beta - 1) * log(m - a)
        return term_1 - term_2
    else
        return -Inf
    end
end

function Distributions.quantile(d::PERT, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    p == 0 && return d.a
    p == 1 && return d.m
    a, b, m = d.a, d.b, d.m
    alpha, beta = _α(d), _β(d)
    value = SpecialFunctions.beta_inc_inv(alpha, beta, p)
    return value[1] * (m - a) + a
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::PERT)
    a, b, m = d.a, d.b, d.m
    alpha, beta = _α(d), _β(d)
    return rand(rng, Distributions.Beta(alpha, beta)) * (m - a) + a
end