"""
    Bradford(a)

A *Bradford* distribution is defined by the following probability density function (PDF):

```math
f(x; a) = \\frac{a}{\\log(1+a) \\cdot (1 + a \\cdot x)}, \\quad 0 < x < 1
```
where:

```julia
Bradford()        # equivalent to Bradford(1)

params(d)        # Get the parameters, i.e. a
```

External links:

* [Bradford on Wikipedia](https://en.wikipedia.org/wiki/Bradford%27s_law)
"""
struct Bradford{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T
    Bradford{T}(a) where {T<:Real} = new{T}(a)
end

# Constructor for integer inputs
Bradford(a::Integer; check_args::Bool=true) = Bradford(float(a); check_args=check_args)

# Main constructor with argument checking
function Bradford(a::Real; check_args::Bool=true)
    @check_args Bradford (a, a > zero(a))
    return Bradford{typeof(a)}(a)
end

Bradford() = Bradford{Float64}(1)
@distr_support Bradford 0.0 1.0

# Parameters
params(d::Bradford) = (d.a,)
@inline partype(d::Bradford{T}) where {T<:Real} = T

Base.eltype(::Type{Bradford{T}}) where {T} = T
# Accessors for individual parameters

# Some statistics 
Statistics.mean(d::Bradford) = (d.a - LogExpFunctions.log1p(d.a))/(d.a * LogExpFunctions.log1p(d.a))
function Statistics.var(d::Bradford) 
    a = d.a
    E_X2 = (d.a^2 / 2 - d.a + LogExpFunctions.log1p(d.a)) / (d.a^2 * LogExpFunctions.log1p(d.a))
    return E_X2 - ((d.a - LogExpFunctions.log1p(d.a))/(d.a * LogExpFunctions.log1p(d.a)))^2
end 

Statistics.median(d::Bradford) = LogExpFunctions.expm1(0.5 * LogExpFunctions.log1p(d.a)) / d.a
StatsBase.entropy(d::Bradford) = log(1+d.a)/2.0 - log(d.a/log(1+d.a))
# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::Bradford, x::Real)
    a = d.a
    isnan(float(x)) && return NaN
    x <= 0 && return 0.0
    x >= 1 && return 1.0

    return LogExpFunctions.log1p(a * x) / LogExpFunctions.log1p(a)
end

function Distributions.pdf(d::Bradford, x::Real)
    a = d.a
    _insupport = insupport(d, x)
    if _insupport
        return a / (LogExpFunctions.log1p(a) * (1 + a*x))
    else
        return 0.0            
    end
end

function Distributions.logpdf(d::Bradford, x::Real)
    a = d.a
    _insupport = insupport(d, x)
    if _insupport
        return log(a) - log(LogExpFunctions.log1p(a)) - log(1 + a*x)
    else
        return -Inf
    end
end

function Distributions.quantile(d::Bradford, p::Real)
    a = d.a
    return LogExpFunctions.expm1(p * LogExpFunctions.log1p(a))/a
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Bradford)
    u = rand(rng)
    return Distributions.quantile(d, u)
end