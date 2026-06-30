"""
    FlorySchulz(a)

A *Flory-Schulz* distribution, commonly used in polymer chemistry to describe the distribution of chain lengths, is defined by the following probability mass function (PMF):

```math
P(X = k) = a^2 k (1 - a)^{k-1}, \\quad k \\in \\{1, 2, 3, \\dots\\}
```
where:

```julia
FlorySchulz()        # equivalent to FlorySchulz(0.5)

params(d)        # Get the parameters, i.e. a
```

External link:

* [Flory Schulz distribution on Wikipedia](https://en.wikipedia.org/wiki/Flory%E2%80%93Schulz_distribution)
"""
struct FlorySchulz{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    a::T
    FlorySchulz{T}(a) where {T<:Real} = new{T}(a)
end

# Constructor for integer inputs
FlorySchulz(a::Integer; check_args::Bool=true) = FlorySchulz(float(a); check_args=check_args)

# Main constructor with argument checking
function FlorySchulz(a::Real; check_args::Bool=true)
    @check_args FlorySchulz (a, 0.0 < a < 1.0)
    return FlorySchulz{typeof(a)}(a)
end

FlorySchulz() = FlorySchulz{Float64}(0.5)
@distr_support FlorySchulz 1 Inf

# Parameters
params(d::FlorySchulz) = (d.a,)
@inline partype(d::FlorySchulz{T}) where {T<:Real} = T

Base.eltype(::Type{FlorySchulz{T}}) where {T} = T
# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::FlorySchulz) = 2/d.a - 1
Statistics.var(d::FlorySchulz) = (2 - 2 * d.a) / d.a^2
StatsBase.mode(d::FlorySchulz) = -(1/log(1 - d.a))
StatsBase.skewness(d::FlorySchulz) = (2 - d.a)/sqrt(2 - 2*d.a)
StatsBase.kurtosis(d::FlorySchulz) = ((d.a - 6)*d.a + 6)/(2 - 2*d.a)
# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::FlorySchulz, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)
    
    return 1 - (1 - a)^x * (1 + a*x)
end

function Distributions.pdf(d::FlorySchulz, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)

    return a^2 * x * (1 - a)^(x-1)
end

function Distributions.logpdf(d::FlorySchulz, x::Real)
    a = d.a
    if !insupport(d, x)
        return -Inf
    end
    x = floor(Int, x)
    return 2 * log(a) + log(x) + (x-1)*log(1 - a) 
end

function Distributions.quantile(d::FlorySchulz, p::Real)
    a = d.a
    if p < 0.0 || p > 1.0
        throw(ArgumentError("p must be in [0, 1]"))
    end

    low = 1
    high = 10

    while 1 - (1 - a)^high * (1 + a * high) < p
        high *= 2
    end

    while low < high
        mid = div(low + high, 2)
        if 1 - (1 - a)^mid * (1 + a * mid) < p
            low = mid + 1
        else
            high = mid
        end
    end

    return low
end

# Sampling
function Distributions.rand(rng::Distributions.AbstractRNG, d::FlorySchulz)
    u = rand(rng)
    return Distributions.quantile(d, u) 
end