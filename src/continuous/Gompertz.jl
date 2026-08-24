"""
    Gompertz(η , b)

A `Gompertz` distribution is commonly used to describe the distribution of adult lifespans and model survival data. The probability density function (PDF) of the Gompertz distribution is given by:

```math
f(x; \\eta, b) = \\eta b \\exp(\\eta + b x - \\eta e^{b x}), \\quad x \\geq 0
```

```julia
Gompertz()        # equivalent to Gompertz(1, 1)
Gompertz(η)        # equivalent to Gompertz(η, 1)

params(d)        # Get the parameters, i.e. (η, b)
```

External links:

* [Gompertz distribution on Wikipedia](https://en.wikipedia.org/wiki/Gompertz_distribution)
"""
struct Gompertz{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    η::T
    b::T
    Gompertz{T}(η, b) where {T<:Real} = new{T}(η, b)
end

Gompertz(η::Integer, b::Integer; check_args::Bool=true) = Gompertz(float(η), float(b); check_args=check_args)

function Gompertz(η::Real, b::Real; check_args::Bool=true)
    @check_args Gompertz (η, η > zero(η)) (b, b > zero(b))
    return Gompertz{typeof(η)}(η, b)
end

Gompertz() = Gompertz{Float64}(1, 1)
Gompertz(η) = Gompertz{Float64}(η, 1.0)
@distr_support Gompertz 0.0 Inf

# parameters

params(d::Gompertz) = (d.η, d.b)
@inline partype(d::Gompertz{T}) where {T<:Real} = T

Base.eltype(::Type{Gompertz{T}}) where {T} = T

#Statistic


Statistics.mean(d::Gompertz) = (1.0/d.b) * exp(d.η) * -SpecialFunctions.expinti(-d.η)

function Statistics.var(d::Gompertz)
    η, b = d.η, d.b

    term_1 = -η * HypergeometricFunctions.pFq((1.0, 1.0, 1.0), (2.0, 2.0, 2.0), -η)
    term_2 = 0.5 * ((π^2 / 6) + (Base.MathConstants.eulergamma + log(η))^2)
    term_3 = (1.0/b) * exp(η) * (-SpecialFunctions.expinti(-η))

    return (2/b^2) * exp(η) * (term_1 + term_2) - term_3^2
end

Statistics.median(d::Gompertz) = (1/d.b) * log(1.0 - 1/d.η *log(0.5))

#### evaluate functions CDF, PDF, logPDF, Quantile, MGF

function Distributions.cdf(d::Gompertz, x::Real)
    η, b = d.η, d.b
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0

    return -expm1(-η * expm1(b * x))
end

function Distributions.ccdf(d::Gompertz, x::Real)
    η, b = d.η, d.b
    isnan(float(x)) && return NaN
    x <= 0 && return 1.0
    isinf(x) && return 0.0

    return exp(-η * expm1(b * x))
end

function Distributions.pdf(d::Gompertz, x::Real)
    η, b = d.η, d.b
    _insupport = insupport(d, x)
    if _insupport
        return b * η * exp(η + b*x - η * exp(b*x))
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Gompertz, x::Real)
    η, b = d.η, d.b
    _insupport = insupport(d, x)
    if _insupport
        return log(b) + log(η) + (η + b*x - η * exp(b*x))
    else
        return -Inf
    end
end

function Distributions.mgf(d::Gompertz, t)
    η, b = d.η, d.b
    return η * exp(η) * SpecialFunctions.expint(t/b, η)
end

function Distributions.quantile(d::Gompertz, p::Real)
    η, b = d.η, d.b
    return (1/b) * log(1 - (1/η)*log(1 - p))
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Gompertz)
    u = rand(rng)
    return Distributions.quantile(d, u)
end