"""
    Benktander_Type2(a, b)

A *Benktander_Type2* distribution, which is used in actuarial science for modeling heavy-tailed distributions. The probability density function (pdf) of the Benktander Type 2 distribution is given by:

```math
f(x; a, b) = \\frac{e^{\\frac{a}{b}(1 - x^b)} \\cdot x^{b-2}}{b} \\left( a \\cdot x^b - b + 1 \\right), \\quad x \\geq 1
```

```julia
Benktander_Type2()        # equivalent to Benktander_Type2(1, 1)
Benktander_Type2(a)        # equivalent to Benktander_Type2(a, 1)

params(d)        # Get the parameters, i.e. (a, b)
```

External links:

* [Benktander Type II distribution on Wikipedia](https://en.wikipedia.org/wiki/Benktander_type_II_distribution)
"""
struct Benktander_Type2{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T
    b::T
    Benktander_Type2{T}(a, b) where {T<:Real} = new{T}(a, b)
end

Benktander_Type2(a::Integer, b::Integer; check_args::Bool=true) = Benktander_Type2(float(a), float(b); check_args=check_args)

function Benktander_Type2(a::Real, b::Real; check_args::Bool=true)
    @check_args Benktander_Type2 (a, a > zero(a)) (b, one(b) >= b > zero(b))
    return Benktander_Type2{typeof(a)}(a, b)
end

Benktander_Type2() = Benktander_Type2{Float64}(1, 1)
Benktander_Type2(a) = Benktander_Type2{Float64}(a, 1.0)
@distr_support Benktander_Type2 1.0 Inf

# parameters

params(d::Benktander_Type2) = (d.a, d.b)
@inline partype(d::Benktander_Type2{T}) where {T<:Real} = T

Base.eltype(::Type{Benktander_Type2{T}}) where {T} = T

#Statistic


Statistics.mean(d::Benktander_Type2) = 1.0 + 1/d.a

function Statistics.var(d::Benktander_Type2)
    a, b = d.a, d.b
    num = -b + 2*a*exp(a/b)*SpecialFunctions.expint(1 - 1/b, a/b)
    dem = a^2 * b
    return num/dem
end

function Statistics.median(d::Benktander_Type2)
    a, b = d.a, d.b
    if b == 1
        return log(2)/a + 1.0
    else
        term_1 = (1-b)/a
        term_2 = (2^(b/(1-b)) * a * exp(a/(1-b)))/(1-b)
        return (term_1 * LambertW.lambertw(term_2))^(1/b)
    end
end

StatsBase.mode(d::Benktander_Type2) = 1.0
#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Benktander_Type2, x::Real)
    a, b = d.a, d.b
    isnan(float(x)) && return NaN
    x <= 1 && return 0.0
    isinf(x) && return 1.0

    return 1 - x^(b - 1) * exp((a / b) * (1 - x^b))
end

function Distributions.pdf(d::Benktander_Type2, x::Real)
    a, b = d.a, d.b
    _insupport = insupport(d, x)
    if _insupport
        return exp((a/b)*(1-x^b)) * x^(b-2) * (a*x^b - b + 1)
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Benktander_Type2, x::Real)
    a, b = d.a, d.b
    _insupport = insupport(d, x)
    if _insupport
        term_1 = (a/b)*(1-x^b)
        term_2 = (b-2)*log(x)
        term_3 = log(a*x^b - b + 1)
        return term_1 + term_2 + term_3
    else
        return -Inf
    end
end

function Distributions.quantile(d::Benktander_Type2, p::Real)
    a, b = d.a, d.b
    if b == 1
        return 1 - log(1-p)/a
    else
        term_1 = (1-b)/a
        term_2 = (a * exp(a/(1-b)) * (1-p)^(-b/(1-b)))/(1-b)
        return (term_1 * LambertW.lambertw(term_2))^(1/b) 
    end
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Benktander_Type2)
    u = rand(rng)
    return Distributions.quantile(d, u)
end