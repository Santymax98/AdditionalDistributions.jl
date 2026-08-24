"""
    Benktander_Type1(a, b)

A *Benktander_Type1* distribution, which is used in actuarial science for modeling heavy-tailed distributions. The probability density function (pdf) of the Benktander Type 1 distribution is given by:

```math
f(x; a, b) = \\left[\\left(1 + \\frac{2b \\log x}{a}\\right)\\left(1 + a + 2b \\log x\\right) - \\frac{2b}{a}\\right] x^{-\\left(2 + a + b \\log x\\right)}, \\quad x \\geq 1
```

```julia
Benktander_Type1()        # equivalent to Benktander_Type1(1, 1)
Benktander_Type1(a)        # equivalent to Benktander_Type1(a, a(a+1)/2)

params(d)        # Get the parameters, i.e. (a, b)
```

External links:

* [Benktander Type I distribution on Wikipedia](https://en.wikipedia.org/wiki/Benktander_type_I_distribution)
"""
struct Benktander_Type1{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T
    b::T
    Benktander_Type1{T}(a, b) where {T<:Real} = new{T}(a, b)
end

Benktander_Type1(a::Integer, b::Integer; check_args::Bool=true) = Benktander_Type1(float(a), float(b); check_args=check_args)

function Benktander_Type1(a::Real, b::Real; check_args::Bool=true)
    @check_args Benktander_Type1 (a, a > zero(a)) (b, a*(a+1)/2 >= b > zero(b))
    return Benktander_Type1{typeof(a)}(a, b)
end

Benktander_Type1() = Benktander_Type1{Float64}(1, 1)
Benktander_Type1(a) = Benktander_Type1{Float64}(a, a*(a+1)/2)
@distr_support Benktander_Type1 1.0 Inf

# parameters

params(d::Benktander_Type1) = (d.a, d.b)
@inline partype(d::Benktander_Type1{T}) where {T<:Real} = T

Base.eltype(::Type{Benktander_Type1{T}}) where {T} = T

#Statistic


Statistics.mean(d::Benktander_Type1) = 1.0 + 1/d.a

function Statistics.var(d::Benktander_Type1)
    a, b = d.a, d.b
    num = -sqrt(b) + a * exp((a-1)^2 / (4*b)) * sqrt(π) * SpecialFunctions.erfc((a-1)/(2*sqrt(b)))
    dem = a^2 * sqrt(b)
    return num/dem
end

#### evaluate functions CDF, PDF, logPDF, quantile

function Distributions.cdf(d::Benktander_Type1, x::Real)
    a, b = d.a, d.b
    isnan(float(x)) && return NaN
    x <= 1 && return 0.0
    isinf(x) && return 1.0

    lx = log(x)
    return 1 - (1 + (2 * b * lx) / a) * exp(-(a + 1 + b * lx) * lx)
end

function Distributions.pdf(d::Benktander_Type1, x::Real)
    a, b = d.a, d.b
    _insupport = insupport(d, x)
    if _insupport
        return (((1 + (2*b*log(x))/a) * (1 + a + 2*b*log(x))) - 2*b/a) * x^(-(2+a+b*log(x)))
    else
        return 0.0
    end
end

function Distributions.logpdf(d::Benktander_Type1, x::Real)
    a, b = d.a, d.b
    _insupport = insupport(d, x)
    if _insupport
        term_1 = log(((1 + (2*b*log(x))/a) * (1 + a + 2*b*log(x))) - 2*b/a)
        term_2 = -(2+a+b*log(x))*log(x)
        return term_1 + term_2
    else
        return -Inf
    end
end

function Distributions.quantile(d::Benktander_Type1, p::Real)
    cdf_func(x) = Distributions.cdf(d, x) - p
    return x = Roots.find_zero(cdf_func, [1.0, 1e6], Roots.Brent())
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Benktander_Type1)
    u = rand(rng)
    return Distributions.quantile(d, u)
end