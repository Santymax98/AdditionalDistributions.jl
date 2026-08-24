"""

    AsymLaplace(m, λ, κ)



An `AsymLaplace` distribution is an asymmetric version of the Laplace

distribution with location parameter `m`, rate parameter `λ`, and asymmetry

parameter `κ`.



The probability density function (PDF) of the AsymLaplace distribution is given by:



```math

f(x; m, λ, κ) = \\frac{λ}{κ + 1/κ}\\begin{cases}
\\exp\\left(\\frac{λ}{κ}(x - m)\\right), & x < m, \\\\

\\exp\\left(-λκ(x - m)\\right), & x \\geq m,

\\end{cases}

\\quad λ > 0, \\quad κ > 0.

```


When `κ = 1`, the distribution reduces to a symmetric Laplace distribution

centered at `m` with rate parameter `λ`.


```julia

AsymLaplace()             # equivalent to AsymLaplace(0, 1, 1)
AsymLaplace(m)            # equivalent to AsymLaplace(m, 1, 1)
AsymLaplace(m, λ)         # equivalent to AsymLaplace(m, λ, 1)

params(d)                 # Get the parameters, i.e. (m, λ, κ)
location(d)               # Get the location parameter m
rate(d)                   # Get the rate parameter λ
asymmetry(d)              # Get the asymmetry parameter κ

```


External links:


* [Asymmetric Laplace distribution on Wikipedia](https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution)
* [Asymmetric Laplace distribution in SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.laplace_asymmetric.html)
* [Kozubowski and Podgórski (2000)](https://doi.org/10.1007/PL00022717)
"""
struct AsymLaplace{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    m::T
    λ::T
    κ::T
    AsymLaplace{T}(m, λ, κ) where {T<:Real} = new{T}(m, λ, κ)
end

function AsymLaplace(m::Real, λ::Real, κ::Real; check_args::Bool=true)
    @check_args AsymLaplace (λ, λ > zero(λ)) (κ, κ > zero(κ))
    m, λ, κ = promote(m, λ, κ)
    return AsymLaplace{typeof(m)}(m, λ, κ) 
end

AsymLaplace(m::Real, λ::Real; check_args::Bool=true) = AsymLaplace(m, λ, 1.0; check_args=check_args)
AsymLaplace(m::Real; check_args::Bool=true) = AsymLaplace(m, 1.0, 1.0; check_args=check_args)
AsymLaplace() = AsymLaplace(0.0, 1.0, 1.0)
AsymLaplace(m::Integer, λ::Integer, κ::Integer; check_args::Bool=true) = AsymLaplace(float(m), float(λ), float(κ); check_args=check_args)
@distr_support AsymLaplace -Inf Inf

params(d::AsymLaplace) = (d.m, d.λ, d.κ)
@inline partype(d::AsymLaplace{T}) where {T<:Real} = T
Base.eltype(::Type{AsymLaplace{T}}) where {T} = T

location(d::AsymLaplace) = d.m
rate(d::AsymLaplace) = d.λ
asymmetry(d::AsymLaplace) = d.κ

Statistics.mean(d::AsymLaplace) = d.m + (1 - d.κ^2) / (d.λ * d.κ)
Statistics.median(d::AsymLaplace) = d.κ < one(d.κ) ? d.m - (1 / (d.λ * d.κ)) * log((1 + d.κ^2) / 2) : d.m + (d.κ / d.λ) * log((1 + d.κ^2) / (2d.κ^2))
Statistics.var(d::AsymLaplace) = (1 + d.κ^4) / (d.λ^2 * d.κ^2)
StatsBase.skewness(d::AsymLaplace) = 2 * (1 - d.κ^6) / ((1 + d.κ^4)^(3/2))
StatsBase.kurtosis(d::AsymLaplace) = 6 * (1 + d.κ^8) / ((1 + d.κ^4)^2)
StatsBase.entropy(d::AsymLaplace) = log(ℯ * (1 + d.κ^2) / (d.λ * d.κ))

function Distributions.cdf(d::AsymLaplace, x::Real)
    m, λ, κ = d.m, d.λ, d.κ
    isnan(float(x)) && return NaN
    isinf(x) && return x > 0 ? 1.0 : 0.0
    if x > m
        return 1 - (1 / (1 + κ^2)) * exp(-λ * κ * (x - m))
    else
        return (κ^2 / (1 + κ^2)) * exp(λ * (x - m) / κ)
    end
end

function Distributions.pdf(d::AsymLaplace, x::Real)
    m, λ, κ = d.m, d.λ, d.κ
    insupport(d, x) || return 0.0
    isinf(x) && return 0.0
    c = λ/(κ + 1/κ)
    if x < m
        return c * exp(λ * (x - m) / κ)
    else
        return c * exp(-λ * κ * (x - m))
    end
end

function Distributions.logpdf(d::AsymLaplace, x::Real)
    m, λ, κ = d.m, d.λ, d.κ

    isnan(float(x)) && return NaN
    insupport(d, x) || return -Inf
    isinf(x) && return -Inf

    logc = log(λ) + log(κ) - log1p(κ^2)

    if x < m
        return logc + λ * (x - m) / κ
    else
        return logc - λ * κ * (x - m)
    end
end

function Distributions.quantile(d::AsymLaplace, q::Real)
    m, λ, κ = d.m, d.λ, d.κ

    (0 <= q <= 1) || throw(DomainError(q, "q must be in [0, 1]"))
    q == 0 && return oftype(float(m), -Inf)
    q == 1 && return oftype(float(m), Inf)

    q0 = κ^2 / (1 + κ^2)

    if q <= q0
        return m + (κ / λ) * log(q * (1 + κ^2) / κ^2)
    else
        return m - (1 / (λ * κ)) * log((1 - q) * (1 + κ^2))
    end
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::AsymLaplace,)
    return quantile(d, rand(rng))
end