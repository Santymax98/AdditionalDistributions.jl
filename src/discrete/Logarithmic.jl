"""
    Logarithmic(a)

A `Logarithmic` distribution, also known as the log-series distribution, is defined by the following probability mass function (PMF):

```math
P(X = k) = -\\frac{1}{\\log(1 - p)} \\cdot \\frac{p^k}{k}, \\quad k \\in \\{1, 2, 3, \\dots\\}
```
where:

```julia
Logarithmic()        # equivalent to Logarithmic(0.5)

params(d)        # Get the parameters, i.e. a
```

External link:

* [Logarithmic distribution on Wikipedia](https://en.wikipedia.org/wiki/Logarithmic_distribution)
"""
struct Logarithmic{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    a::T
    Logarithmic{T}(a) where {T<:Real} = new{T}(a)
end

# Constructor for integer inputs
Logarithmic(a::Integer; check_args::Bool=true) = Logarithmic(float(a); check_args=check_args)

# Main constructor with argument checking
function Logarithmic(a::Real; check_args::Bool=true)
    @check_args Logarithmic (a, 0.0 < a < 1.0)
    return Logarithmic{typeof(a)}(a)
end

Logarithmic() = Logarithmic{Float64}(0.5)
@distr_support Logarithmic 1 Inf

# Parameters
params(d::Logarithmic) = (d.a,)
@inline partype(d::Logarithmic{T}) where {T<:Real} = T

Base.eltype(::Type{Logarithmic{T}}) where {T} = T
# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::Logarithmic) = (-1/log(1-d.a)) * (d.a/(1-d.a))
Statistics.var(d::Logarithmic) = -(d.a^2 + d.a * log(1-d.a))/((1-d.a)^2 * (log(1-d.a))^2)
StatsBase.mode(d::Logarithmic) = 1

# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::Logarithmic, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)
    ks = 1:x
    s = sum(a .^ ks ./ ks)
    return -s / log(1 - a)
end

function Distributions.pdf(d::Logarithmic, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)

    return -a^x/(x*log(1-a))
end

function Distributions.logpdf(d::Logarithmic, x::Real)
    a = d.a
    if !insupport(d, x)
        return -Inf
    end
    x = floor(Int, x)
    return log(-1/(log(1-a)) * (a^x / x))
end

function Distributions.quantile(d::Logarithmic, p::Real)
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
function Distributions.rand(rng::Distributions.AbstractRNG, d::Logarithmic{T}) where T
    a = d.a
    h = log(1 - a)
    u2 = rand(rng, Distributions.Uniform(0, 1))
    x = 1
    if u2 > a
        return x
    else
        u1 = rand(rng, Distributions.Uniform(0, 1))
        q = 1 - exp(u1*h)
        if u2 < q^2
            return Int(trunc(1 + log(u2) / log(q)))
        else
            if u2 > q
                return 1
            else
                return 2
            end
        end
    end
end