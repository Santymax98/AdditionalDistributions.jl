"""
    HalfCauchy(σ)

A `HalfCauchy` distribution is the distribution of the absolute value of a
zero-centered Cauchy random variable. Equivalently, if ``Z \\sim Cauchy(0, σ)``,
then ``X = |Z|`` follows a `HalfCauchy(σ)` distribution.

The probability density function (PDF) of the HalfCauchy distribution is given by:

```math 
f(x; σ) = \\frac{2}{π σ}\\frac{1}{1 + (x/σ)^2}, \\quad x \\geq 0, \\quad σ > 0.
```

```julia
HalfCauchy()         # equivalent to HalfCauchy(1)
HalfCauchy(σ)

params(d)            # Get the parameters, i.e. (σ,)
```

External links:

* [Half-Cauchy distribution on Wikipedia](https://distribution-explorer.github.io/continuous/halfcauchy.html)
"""
struct HalfCauchy{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    σ::T
    HalfCauchy{T}(σ) where {T<:Real} = new{T}(σ)
end

HalfCauchy(σ::Integer; check_args::Bool=true) = HalfCauchy(float(σ); check_args=check_args)

function HalfCauchy(σ::Real; check_args::Bool=true)
    @check_args HalfCauchy (σ, σ > zero(σ))
    return HalfCauchy{typeof(σ)}(σ)
end

HalfCauchy() = HalfCauchy(1.0)
@distr_support HalfCauchy 0.0 Inf

params(d::HalfCauchy) = (d.σ,)
@inline partype(d::HalfCauchy{T}) where {T<:Real} = T
Base.eltype(::Type{HalfCauchy{T}}) where {T} = T

scale(d::HalfCauchy) = d.σ

function Distributions.cdf(d::HalfCauchy, x::Real)
    σ = d.σ
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    isinf(x) && return 1.0
    return (2 / π) * atan(x / σ)
end

function Distributions.pdf(d::HalfCauchy, x::Real)
    σ = d.σ
    insupport(d, x) || return 0.0
    isinf(x) && return 0.0
    return (2 / (π * σ)) * (1 / (1 + (x / σ)^2))
end

function Distributions.logpdf(d::HalfCauchy, x::Real)
    σ = d.σ
    insupport(d, x) || return -Inf
    isinf(x) && return -Inf
    return log(2 / (π * σ)) - log(1 + (x / σ)^2)
end

function Distributions.quantile(d::HalfCauchy, q::Real)
    σ = d.σ
    (0 <= q <= 1) || throw(DomainError(q, "q must be in [0, 1]"))
    q == 0 && return 0.0
    q == 1 && return Inf
    return σ * tan((π / 2) * q)
    
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::HalfCauchy)
    σ = d.σ
    return abs(rand(rng, Distributions.Cauchy(0, σ)))
end