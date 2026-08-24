"""
    HalfTDist(ν, σ)

A `HalfTDist` distribution is the distribution of the absolute value of a
zero-centered Student's t random variable with `ν` degrees of freedom and scale
parameter `σ`. Equivalently, if ``T \\sim t_ν``, then
``X = σ |T|`` follows a `HalfTDist(ν, σ)` distribution.

The probability density function (PDF) of the HalfTDist distribution is given by:

```math
f(x; ν, σ) =
\\frac{2}{σ}
\\frac{Γ\\left(\\frac{ν+1}{2}\\right)}
{\\sqrt{νπ}\\, Γ\\left(\\frac{ν}{2}\\right)}
\\left(1 + \\frac{x^2}{νσ^2}\\right)^{-\\frac{ν+1}{2}},
\\quad x \\geq 0, \\quad ν > 0, \\quad σ > 0.
```

```julia
HalfTDist()          # equivalent to HalfTDist(1, 1)
HalfTDist(ν)         # equivalent to HalfTDist(ν, 1)

params(d)            # Get the parameters, i.e. (ν, σ)
```

External links:

* [Half-t distribution](https://distribution-explorer.github.io/continuous/halft.html)
* [Student's t-distribution on Wikipedia](https://en.wikipedia.org/wiki/Student%27s_t-distribution)
"""
struct HalfTDist{T<:Real} <: Distributions.ContinuousUnivariateDistribution
  ν::T
  σ::T
  HalfTDist{T}(ν, σ) where {T<:Real} = new{T}(ν, σ)
end

function HalfTDist(ν::Real, σ::Real; check_args::Bool=true)
    ν, σ = promote(float(ν), float(σ))
    @check_args HalfTDist (ν, ν > zero(ν)) (σ, σ > zero(σ))
    return HalfTDist{typeof(ν)}(ν, σ)
end

HalfTDist(ν::Real; check_args::Bool=true) = HalfTDist(ν, one(float(ν)); check_args=check_args)
HalfTDist(ν::Integer, σ::Integer; check_args::Bool=true) = HalfTDist(float(ν), float(σ); check_args=check_args)
HalfTDist() = HalfTDist(1.0, 1.0)

@distr_support HalfTDist 0 Inf
params(d::HalfTDist) = (d.ν, d.σ)
@inline partype(d::HalfTDist{T}) where {T<:Real} = T
Base.eltype(::Type{HalfTDist{T}}) where {T} = T

dof(d::HalfTDist) = d.ν
scale(d::HalfTDist) = d.σ

function Statistics.mean(d::HalfTDist)
    ν, σ = d.ν, d.σ
    ν <= 1 && return Inf

    logμ =
        log(σ) +
        log(2) +
        0.5 * log(ν) +
        SpecialFunctions.loggamma((ν + 1) / 2) -
        log(ν - 1) -
        0.5 * log(π) -
        SpecialFunctions.loggamma(ν / 2)

    return exp(logμ)
end

function Statistics.var(d::HalfTDist)
    ν, σ = d.ν, d.σ
    ν <= 2 && return Inf

    μ = mean(d)
    return σ^2 * ν / (ν - 2) - μ^2
end
StatsBase.mode(d::HalfTDist) = zero(d.σ)
Statistics.median(d::HalfTDist) = d.σ * quantile(Distributions.TDist(d.ν), 0.75)

function Distributions.cdf(d::HalfTDist, x::Real)
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0

    ν, σ = d.ν, d.σ
    z = x / σ
    td = Distributions.TDist(ν)

    return 1 - 2 * ccdf(td, z)
end

function Distributions.pdf(d::HalfTDist, x::Real)
    insupport(d, x) || return zero(float(x))
    isinf(x) && return zero(float(x))

    ν, σ = d.ν, d.σ
    z = x / σ

    return 2 * pdf(Distributions.TDist(ν), z) / σ
end

function Distributions.logpdf(d::HalfTDist, x::Real)
    insupport(d, x) || return -Inf
    isinf(x) && return -Inf

    ν, σ = d.ν, d.σ
    z = x / σ

    return log(2) - log(σ) + Distributions.logpdf(Distributions.TDist(ν), z)
end

function Distributions.quantile(d::HalfTDist, q::Real)
    ν, σ = d.ν, d.σ
    (0 <= q <= 1) || throw(DomainError(q, "q must be in [0, 1]"))
    q == 0 && return zero(σ)
    q == 1 && return oftype(float(σ), Inf)
    td = Distributions.TDist(ν)
    if q < 0.5
        return σ * quantile(td, (q + 1) / 2)
    else
        return σ * cquantile(td, (1 - q) / 2)
    end
end

function Distributions.rand(rng::Random.AbstractRNG, d::HalfTDist)
    return d.σ * abs(rand(rng, Distributions.TDist(d.ν)))
end

function Distributions.ccdf(d::HalfTDist, x::Real)
    isnan(float(x)) && return NaN
    x <= 0 && return 1.0
    isinf(x) && return 0.0

    ν, σ = d.ν, d.σ
    z = x / σ
    td = Distributions.TDist(ν)

    return 2 * ccdf(td, z)
end