"""
    Nakagami(m, Ω)

A *Nakagami* distribution, often used to model physical phenomena in areas such as medical ultrasound imaging and communications, is defined by the following probability density function (PDF):

```math
f(x; m, \\Omega) = \\frac{2m^m}{\\Gamma(m) \\Omega^m} x^{2m-1} \\exp\\left(-\\frac{m}{\\Omega} x^2\\right), \\quad x \\geq 0
```
where:

- m is a shape parameter
- ``\\Omega`` is a spread parameter

```julia
Nakagami()        # equivalent to Nakagami(0.5, 1)
Nakagami(m)        # equivalent to Nakagami(m, 1)

params(d)        # Get the parameters, i.e. (m, Ω)
```

External links:

* [Nakagami distribution on Wikipedia](https://en.wikipedia.org/wiki/Nakagami_distribution)
"""
struct Nakagami{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    m::T
    Ω::T
    Nakagami{T}(m, Ω) where {T<:Real} = new{T}(m, Ω)
end

Nakagami(m::Integer, Ω::Integer; check_args::Bool=true) = Nakagami(float(m), float(Ω); check_args=check_args)

function Nakagami(m::Real, Ω::Real; check_args::Bool=true)
    @check_args Nakagami (m, m >= 0.5) (Ω, Ω > zero(Ω))
    return Nakagami{typeof(m)}(m, Ω)
end

Nakagami() = Nakagami{Float64}(0.5, 1)
Nakagami(m) = Nakagami{Float64}(m, 1)
@distr_support Nakagami 0.0 Inf

# parameters
params(d::Nakagami) = (d.m, d.Ω)
@inline partype(d::Nakagami{T}) where {T<:Real} = T

Base.eltype(::Type{Nakagami{T}}) where {T} = T

shape(d::Nakagami) = d.m
spread(d::Nakagami) = d.Ω

#Statistic

Statistics.mean(d::Nakagami) = (SpecialFunctions.gamma(d.m + 0.5)/SpecialFunctions.gamma(d.m)) * (d.Ω/d.m)^(1/2)
Statistics.var(d::Nakagami) = d.Ω * (1 - (1/d.m) * ((SpecialFunctions.gamma(d.m + 0.5))/SpecialFunctions.gamma(d.m))^2)
Statistics.median(d::Nakagami) = Distributions.quantile(d, 0.5)
StatsBase.mode(d::Nakagami) = ((2*d.m - 1)*d.Ω / (2*d.m))^(1/2)
#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Nakagami, x::Real)
    m, Ω = d.m, d.Ω
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0

    return SpecialFunctions.gamma_inc(m, (m / Ω) * x^2)[1]
end

function Distributions.pdf(d::Nakagami, x::Real)
    m, Ω = d.m, d.Ω
    _insupport = insupport(d, x)
    if _insupport
        term1 = (2*m^m / (SpecialFunctions.gamma(m) * Ω^m)) * x^(2*m - 1)
        term2 = exp(-(m/Ω) * x^2)
        return term1 * term2
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Nakagami, x::Real)
    m, Ω = d.m, d.Ω
    _insupport = insupport(d, x)
    if _insupport
        term_1 = log(2 * m^m) - SpecialFunctions.loggamma(m) - m * log(Ω) + (2 * m - 1) * log(x)
        term_2 = -(m / Ω) * x^2
        return term_1 + term_2
    else
        return -Inf
    end
end

function Distributions.quantile(d::Nakagami, p::Real)
    m, Ω = d.m, d.Ω
    z =  SpecialFunctions.gamma_inc_inv(m, p, 1 - p)
    return sqrt(Ω / m * z)
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Nakagami)
    m, Ω = d.m, d.Ω #Nakagami parameters
    k, θ = m, Ω/m #Gamma parameters
    Y = rand(rng, Distributions.Gamma(k, θ))
    return sqrt(Y)
end