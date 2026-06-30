"""
    Maxwell(a)

A *Maxwell-Boltzmann* distribution, often used in physics to describe the distribution of particle speeds in idealized gases, is defined by the following probability density function (PDF):

```math
f(x; a) = \\sqrt{\\frac{2}{\\pi}} \\frac{x^2}{a^3} \\exp\\left(-\\frac{x^2}{2a^2}\\right), \\quad x > 0
```

```julia
Maxwell()        # equivalent to Maxwell(1)

params(d)        # Get the parameters, i.e. a
```

External links:

* [Maxwell Boltzmann distribution on Wikipedia](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution)
"""
struct Maxwell{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T
    Maxwell{T}(a) where {T<:Real} = new{T}(a)
end

# Constructor for integer inputs
Maxwell(a::Integer; check_args::Bool=true) = Maxwell(float(a); check_args=check_args)

# Main constructor with argument checking
function Maxwell(a::Real; check_args::Bool=true)
    @check_args Maxwell (a, a > zero(a))
    return Maxwell{typeof(a)}(a)
end

Maxwell() = Maxwell{Float64}(1)
@distr_support Maxwell 0.0 Inf

# Parameters
params(d::Maxwell) = (d.a,)
@inline partype(d::Maxwell{T}) where {T<:Real} = T

Base.eltype(::Type{Maxwell{T}}) where {T} = T
# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::Maxwell) = 2 * d.a * sqrt(2/π)
Statistics.var(d::Maxwell) = (d.a^2 * (3 * π - 8)/π)
Statistics.median(d::Maxwell) = d.a * sqrt(2 * SpecialFunctions.gamma_inc_inv(3.0/2.0, 0.5, 0.5)[1])
# Mode
StatsBase.mode(d::Maxwell) = d.a * sqrt(2)
StatsBase.skewness(d::Maxwell) = (2 * sqrt(2) * (16 - 5*π))/(3*π - 8)^(3/2)
StatsBase.kurtosis(d::Maxwell) = (-192 + π*(16 + 15*π))/(3*π - 8)^2
StatsBase.entropy(d::Maxwell) = log(d.a * sqrt(2*π)) + Base.MathConstants.eulergamma - 0.5
# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::Maxwell, x::Real)
    a = d.a
    _insupport = insupport(d, x)
    if _insupport
        return SpecialFunctions.erf(x/(sqrt(2)*a)) - sqrt(2/π) * (x/a) * exp(-(x^2/(2*a^2)))
    else
        return 0.0            
    end
end

function Distributions.pdf(d::Maxwell, x::Real)
    a = d.a
    _insupport = insupport(d, x)
    if _insupport
        return sqrt(2/π) * (x^2 / a^3) * exp(-(x^2/(2*a^2)))
    else
        return 0.0            
    end
end

function Distributions.logpdf(d::Maxwell, x::Real)
    a = d.a
    _insupport = insupport(d, x)
    if _insupport
        return 0.5 * log(2/π) + 2*log(x) - 3*log(a) - (x^2/(2*a^2))
    else
        return -Inf
    end
end

function Distributions.quantile(d::Maxwell, p::Real)
    a = d.a
    return a * sqrt(2 * SpecialFunctions.gamma_inc_inv(3.0/2.0, p, 1-p)[1])
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Maxwell)
    u = rand(rng)
    return Distributions.quantile(d, u)
end