"""
    Yule(a)

The *Yule* distribution is a discrete probability distribution defined by the following probability mass function (PMF):

```math
P(X = k) = a * B(k, a + 1)
```
where: 

- ``B(k, a + 1)`` is a beta function

```julia
Yule()        # equivalent to Yule(1)

params(d)        # Get the parameters, i.e. a
```

External link
* [Yule distribution on Wikipedia](https://en.wikipedia.org/wiki/Yule%E2%80%93Simon_distribution)
"""
struct Yule{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    a::T
    Yule{T}(a) where {T<:Real} = new{T}(a)
end

# Constructor for integer inputs
Yule(a::Integer; check_args::Bool=true) = Yule(float(a); check_args=check_args)

# Main constructor with argument checking
function Yule(a::Real; check_args::Bool=true)
    @check_args Yule (a, zero(a) < a)
    return Yule{typeof(a)}(a)
end

Yule() = Yule{Float64}(1)
@distr_support Yule 1 Inf

# Parameters
params(d::Yule) = (d.a,)
@inline partype(d::Yule{T}) where {T<:Real} = T

Base.eltype(::Type{Yule{T}}) where {T} = T
# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::Yule) = d.a > 1 ? d.a/(d.a - 1) : NaN
Statistics.var(d::Yule) = d.a > 2 ? d.a^2 /((d.a-1)^2 *(d.a-2)) : NaN
StatsBase.mode(d::Yule) = 1

# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::Yule, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)

    return 1 - x * SpecialFunctions.beta(x, a+1)
end

function Distributions.pdf(d::Yule, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)

    return a * SpecialFunctions.beta(x, a+1)
end

function Distributions.logpdf(d::Yule, x::Real)
    a = d.a
    if !insupport(d, x)
        return -Inf
    end
    x = floor(Int, x)

    return log(a) + SpecialFunctions.logbeta(x, a+1)
end

function Distributions.quantile(d::Yule, p::Real)
    a = d.a
    if p < 0.0 || p > 1.0
        throw(ArgumentError("p must be in [0, 1]"))
    end

    low = 1
    high = 10

    while Distributions.cdf(d, high) < p
        high *= 2
    end

    while low < high
        mid = div(low + high, 2)
        if Distributions.cdf(d, mid) < p
            low = mid + 1
        else
            high = mid
        end
    end

    return low
end

# Sampling
function Distributions.rand(rng::Distributions.AbstractRNG, d::Yule{T}) where T
    a = d.a
    W = rand(rng, Distributions.Exponential(1/a))
    return 1 + rand(rng, Distributions.Geometric(exp(-W)))
end