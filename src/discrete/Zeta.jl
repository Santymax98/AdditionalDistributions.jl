"""
    Zeta(s)

The *Zeta distribution* is a discrete probability distribution defined by the following probability mass function (PMF):

```math
P(X = k) = \\frac{1}{\\zeta(s)} \\cdot \\frac{1}{k^s}
```
Where ``\\zeta`` is a Riemann Function

```julia
Zeta()      # equivalent to Zeta(1)

params(d)   # Get the parameters, i.e. s
```
External link:

*[Zeta Distribution on Wikipedia](https://en.wikipedia.org/wiki/Zeta_distribution)
"""
struct Zeta{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    s::T
    Zeta{T}(s) where {T<:Real} = new{T}(s)
end

Zeta(s::Integer; check_args::Bool=true) = Zeta(float(s); check_args=check_args)
function Zeta(s::Real; check_args::Bool=true)
    @check_args Zeta (s, one(s) < s)
    return Zeta{typeof(s)}(s)
end

Zeta() = Zeta(2.0)

@distr_support Zeta 1 Inf

# parameters
params(d::Zeta) = (d.s,)
# statistics 
Statistics.mean(d::Zeta) = d.s > 2 ? SpecialFunctions.zeta(d.s - 1)/SpecialFunctions.zeta(d.s) : NaN
Statistics.var(d::Zeta) = d.s > 3 ? (SpecialFunctions.zeta(d.s) * SpecialFunctions.zeta(d.s - 2) - SpecialFunctions.zeta(d.s - 1)^2)/SpecialFunctions.zeta(d.s)^2 : NaN
#Statistics.median(d::Zeta) = 
StatsBase.mode(d::Zeta) = 1 

#evaluate functions CDF, PDF, logPDF, Quantil
function Distributions.cdf(d::Zeta, x::Real)
    isnan(float(x)) && return NaN
    x == Inf && return 1.0
    !isfinite(x) && return 0.0
    x < 1 && return 0.0

    k = floor(Int, x)
    num = harmonic(k, d.s)
    den = SpecialFunctions.zeta(d.s)
    return num / den
end

@inline _isposinteger_zeta(x::Real) = isfinite(x) && x >= 1 && x == floor(x)

function Distributions.pdf(d::Zeta, x::Real)
    _isposinteger_zeta(x) || return 0.0
    k = floor(Int, x)
    return inv(k^d.s * SpecialFunctions.zeta(d.s))
end

function Distributions.logpdf(d::Zeta, x::Real)
    _isposinteger_zeta(x) || return -Inf
    k = floor(Int, x)
    return -d.s * log(k) - log(SpecialFunctions.zeta(d.s))
end

function Distributions.quantile(d::Zeta, p::Real)
    s = d.s
    if p < 0 || p > 1
        throw(DomainError(p, "p must be in [0, 1]"))
    end
    if p == 0
        return 1
    elseif p == 1
        return Inf
    end

    left, right = 1, 2
    while Distributions.cdf(d, right) < p
        left = right
        right *= 2
    end

    while left < right
        mid = (left + right) ÷ 2
        if Distributions.cdf(d, mid) < p
            left = mid + 1
        else
            right = mid
        end
    end

    return left
end

# sampling
function Distributions.rand(rng::Distributions.AbstractRNG, d::Zeta)
    u = rand(rng)
    return Distributions.quantile(d, u)
end
