"""
    Argus(χ, c)

A *Argus* distribution, which is used in particle physics to describe the invariant mass of a decayed particle candidate. The probability density function (pdf) of the ARGUS distribution is given by:

```math
f(x; χ, c) = \\frac{\\xi^3}{\\sqrt{2\\pi}\\Phi(\\xi)} \\cdot \\frac{x}{c^2} \\sqrt{1 - \\frac{x^2}{c^2}} \\exp\\left(-\\frac{1}{2}\\xi^2\\left(1 - \\frac{x^2}{c^2}\\right)\\right)
```

```julia
Argus()        # equivalent to Argus(1, 1)

params(d)        # Get the parameters, i.e. (χ, c)
```

External links:

* [Argus distribution on Wikipedia](https://en.wikipedia.org/wiki/ARGUS_distribution)
"""
struct Argus{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    χ::T
    c::T
    Argus{T}(χ, c) where {T<:Real} = new{T}(χ, c)
end

Argus(χ::Integer, c::Integer; check_args::Bool=true) = Argus(float(χ), float(c); check_args=check_args)

function Argus(χ::Real, c::Real; check_args::Bool=true)
    @check_args Argus (χ, χ > zero(χ)) (c, c > zero(c))
    return Argus{typeof(χ)}(χ, c)
end

Argus() = Argus{Float64}(1, 1)
@distr_support Argus 0.0 d.c

#aux 
function Ψ(x::Real)
    Distributions.cdf(normal_dist, x) - x * Distributions.pdf(normal_dist, x) - 0.5
end

# parameters

params(d::Argus) = (d.χ, d.c)
@inline partype(d::Argus{T}) where {T<:Real} = T

Base.eltype(::Type{Argus{T}}) where {T} = T

#Statistic


Statistics.mean(d::Argus) = d.c * sqrt(π / 8) *  (d.χ * exp(-d.χ^2 / 4) * SpecialFunctions.besseli(1, d.χ^2 / 4)) / Ψ(d.χ)

function Statistics.var(d::Argus)
    χ, c = d.χ, d.c
    mean_d = Statistics.mean(d)
    
    E_X2 = d.c^2 * (1 - (3/d.χ^2) + ((d.χ * Distributions.pdf(normal_dist, d.χ) / Ψ(d.χ))))
    var_d = E_X2 - mean_d^2
    return var_d
end

StatsBase.mode(d::Argus) = (d.c / (sqrt(2) * d.χ)) * sqrt((d.χ^2 - 2) + sqrt(d.χ^4 + 4))

#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Argus, x::Real)
    χ, c = d.χ, d.c
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    x >= c && return 1.0

    term = χ * sqrt(1 - x^2 / c^2)
    return 1 - Ψ(term) / Ψ(χ)
end

function Distributions.pdf(d::Argus, x::Real)
    χ, c = d.χ, d.c
    _insupport = insupport(d, x)
    if _insupport
        term_1 = χ^3 / (sqrt(2 * π) * Ψ(χ))
        term_2 = (x/c^2) * sqrt(1 - (x^2 / c^2))
        term_3 = exp(-0.5 * χ^2 * (1- (x^2 / c^2)))
        return term_1 * term_2 * term_3
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Argus, x::Real)
    χ, c = d.χ, d.c
    _insupport = insupport(d, x)
    if _insupport
        term_1 = log(χ^3 / (sqrt(2 * π) * Ψ(χ)))
        term_2 = log(x/c^2) + 0.5 * log(1 - (x^2 / c^2))
        term_3 = -0.5 * χ^2 * (1- (x^2 / c^2))
        return term_1 + term_2 + term_3
    else
        return -Inf
    end
end

function Distributions.quantile(d::Argus, p::Real)
    c = d.c

    cdf_func(x) = Distributions.cdf(d, x) - p
    x = Roots.find_zero(cdf_func, [0.0, c] , Roots.Brent())
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Argus)
    u = rand(rng)
    return Distributions.quantile(d, u)
end