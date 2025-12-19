
"""
    IrwinHall(n, a, b)

``IrwinHall(n, a, b)`` distributed variable is a sum of ``n`` independent ``Uniform(a, b)`` variables. Its probability density function (PDF) is given by:

```math
f(x; n, a, b) = \\frac{1}{(b-a)(n-1)!}\\sum_{k=0}^{\\lfloor y\\rfloor}(-1)^k{n \\choose k} (y-k)^{n-1}, \\quad y = \\frac{x - a}{b - a}.
```

Related Bates distribution is ``IrwinHall(n, a/n, b/n)``.

```julia
IrwinHall(n)       # equivalent to IrwinHall(n, 0, 1)

params(d)        # Get the parameters, i.e. (n, a, b)
```

External links:

* [Irwin-Hall distribution on Wikipedia](https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution)
* [Bates distribution on Wikipedia](https://en.wikipedia.org/wiki/Bates_distribution)
"""
struct IrwinHall{S<:Integer, T<:Real} <: Distributions.ContinuousUnivariateDistribution
    n::S
    a::T
    b::T
    IrwinHall{S,T}(n, a, b) where {S<:Integer,T<:Real} = new{S,T}(n, a, b)
end

function IrwinHall(n::Integer, a::T, b::T; check_args::Bool=true) where {T<:Real}
    @check_args IrwinHall (n, n > zero(n)) (b, b > a) 
    return IrwinHall{typeof(n),typeof(a)}(n, a, b)
end

IrwinHall(n::Integer, a::Real, b::Real; check_args::Bool=true) = IrwinHall(n, promote(a, b)...; check_args=check_args)
IrwinHall(n::Integer, a::Integer, b::Integer; check_args::Bool=true) = IrwinHall(n, float(a), float(b); check_args=check_args)

IrwinHall(n::Integer) = IrwinHall(n, zero(n), one(n))

@distr_support IrwinHall (d.n*d.a) (d.n*d.b)

# parameters

params(d::IrwinHall) = (d.n, d.a, d.b)

shape(d::IrwinHall) = d.n
location(d::IrwinHall) = d.n*d.a
scale(d::IrwinHall) = d.b - d.a

# moments, mode, median

Statistics.mean(d::IrwinHall) = d.n*(d.a + d.b)/2
Statistics.var(d::IrwinHall) = d.n*(d.b - d.a)^2 / 12
StatsBase.skewness(d::IrwinHall) = zero(d.a)  
StatsBase.kurtosis(d::IrwinHall) = -6 / (5*d.n)

StatsBase.median(d::IrwinHall) = d.n*(d.a + d.b)/2
StatsBase.mode(d::IrwinHall) = d.n*(d.a + d.b)/2

# pdf, logpdf, cdf, quantile, cf, mgf, cgf

function Distributions.pdf(d::IrwinHall, x::Real)
    n, a, b = d.n, d.a, d.b
    insupport(d, x) || return zero(a)

    n == one(n) && return 1/(b - a) # uniform distribution

    y = (x - n*a) / (b - a)
    y = ifelse(y > n/2, n - y, y) # using symmetry we make calculation faster for large x
    m = floor(Int, y)

    c = one(y)/((b - a) * factorial(n - 1) )
    S = c * y^(n - 1)

    for k in 1:m
        c *= -(n + 1 - k)/k
        S += c * (y - k)^(n - 1)
    end

    return S
end

function Distributions.logpdf(d::IrwinHall, x::Real)
    n, a, b = d.n, d.a, d.b
    insupport(d, x) || return typemin(a)

    y = (x - n*a) / (b - a)
    y <= 1 && return (n-1)*log(y) - SpecialFunctions.logfactorial(n-1) - log(b - a) # left edge
    y >= (n-1) && return (n-1)*log(n - y) - SpecialFunctions.logfactorial(n-1) - log(b - a) # right edge
    return log(pdf(d, x)) # middle
end

function Distributions.cdf(d::IrwinHall, x::Real)
    n, a, b = d.n, d.a, d.b
    x <= minimum(d) && return zero(a)
    x >= maximum(d) && return one(a)

    n == one(n) && return (x - a)/(b - a) # uniform distribution

    y = (x - n*a) / (b - a)
    flag = y > n/2
    y = ifelse(flag, n - y, y) # using symmetry we make calculation faster for large x
    m = floor(Int, y)

    c = one(y)/(n * factorial(n - 1) )
    S = c * y^n

    for k in 1:m
        c *= -(n + 1 - k)/k
        S += c * (y - k)^n
    end

    return flag ? 1-S : S
end

function Distributions.quantile(d::IrwinHall, p::Real)
    n, a, b = d.n, d.a, d.b
    A, B = n*a, n*b
    cdf_func(x) = cdf(d, x) - p
    pdf_func(x) = pdf(d, x)
    initial_point = A + p * (B - A) 

    return Roots.find_zero((cdf_func, pdf_func), initial_point, Roots.Newton())
end

Distributions.cf(d::IrwinHall, t::Real) = cf(Uniform(d.a, d.b), t) ^ d.n
Distributions.mgf(d::IrwinHall, t::Real) = mgf(Uniform(d.a, d.b), t) ^ d.n
Distributions.cgf(d::IrwinHall, t::Real) = d.n * cgf(Uniform(d.a, d.b), t)

# affine transformations, convolution, conversion

Base.:+(d::IrwinHall, x::Real) = IrwinHall(d.n, d.a + x/d.n, d.b + x/d.n)
Base.:*(c::Real, d::IrwinHall) = IrwinHall(d.n, minmax(c*d.a, c*d.b)...)

function Distributions.convolve(d1::IrwinHall, d2::IrwinHall)
    d1.a ≈ d2.a || throw(ArgumentError("$(d1.a) !≈ $(d2.a): a parameters must be approximately equal"))
    d1.b ≈ d2.b || throw(ArgumentError("$(d1.b) !≈ $(d2.b): b parameters must be approximately equal"))

    return IrwinHall(d1.n + d2.n, d1.a, d1.b)
end

Distributions.convert(::Type{IrwinHall}, d::Uniform) = IrwinHall(1, d.a, d.b)

# sampling

Distributions.rand(rng::Distributions.AbstractRNG, d::IrwinHall) = (d.b - d.a)*sum(rand(rng) for _ in 1:d.n) + d.n*d.a
