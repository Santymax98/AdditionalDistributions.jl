"""
    Borel(a)

A `Borel` distribution is a discrete probability distribution often used in branching processes and queueing theory. The probability mass function (PMF) of the Borel distribution is given by:

```math
P(X = k) = \\frac{e^{-a k} (a k)^{k-1}}{k!}, \\quad k \\in \\{1, 2, 3, \\dots\\}
```

```julia
Borel()        # equivalent to Borel(0)

params(d)        # Get the parameters, i.e. a
```

External link:

* [Borel distribution on Wikipedia](https://en.wikipedia.org/wiki/Borel_distribution)
"""
struct Borel{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    a::T
    Borel{T}(a) where {T<:Real} = new{T}(a)
end

# Constructor for integer inputs
Borel(a::Integer; check_args::Bool=true) = Borel(float(a); check_args=check_args)

# Main constructor with argument checking
function Borel(a::Real; check_args::Bool=true)
    @check_args Borel (a, 0.0 <= a <= 1.0)
    return Borel{typeof(a)}(a)
end

Borel() = Borel{Float64}(0)
@distr_support Borel 1 Inf

# Parameters
params(d::Borel) = (d.a,)
@inline partype(d::Borel{T}) where {T<:Real} = T

Base.eltype(::Type{Borel{T}}) where {T} = T
# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::Borel) = 1/(1- d.a)
Statistics.var(d::Borel) = d.a/(1-d.a)^3

# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::Borel, x::Real)
    a = d.a
    isnan(float(x)) && return NaN
    x < 1 && return 0.0
    x == Inf && return 1.0
    !isfinite(x) && return 0.0

    x = floor(Int, x)
    cdf_value = 0.0
    for k in 1:x
        log_pdf = -a * k + (k-1)*log(a * k) - SpecialFunctions.loggamma(k + 1)
        cdf_value +=  exp(log_pdf)
    end
    return cdf_value
end

function Distributions.pdf(d::Borel, x::Real)
    a = d.a
    if !insupport(d, x)
        return 0.0
    end
    x = floor(Int, x)
    pdf_value = -a * x + (x-1)*log(a * x) - SpecialFunctions.loggamma(x + 1)
    return exp(pdf_value)
end

function Distributions.logpdf(d::Borel, x::Real)
    a = d.a
    if !insupport(d, x)
        return -Inf
    end
    x = floor(Int, x)
    return -a * x + (x-1) * log(a * x) - SpecialFunctions.logfactorial(x)
end

function Distributions.quantile(d::Borel, p::Real)
    a = d.a
    if p < 0.0 || p > 1.0
        throw(ArgumentError("p must be in [0,1]"))
    end

    # binary search
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

# Sampling function using rejection sampling with Poisson transformation
function Distributions.rand(rng::Distributions.AbstractRNG, d::Borel)
    u = rand(rng)
    return Distributions.quantile(d, u) 
end