"""
    Bhattacharjee(a, b, σ)

A *Bhattacharjee* distribution is a continuous univariate distribution where the mean follows a rectangular (uniform) distribution and the conditional distribution given the mean is normal. Specifically, if the mean θ follows a uniform distribution between `a` and `b`, and the conditional distribution of `X` given `θ` is normal with mean `θ` and standard deviation `σ`, then `X` follows a Bhattacharjee distribution.

```math
f(x) = \\frac{1}{b - a} \\left[ \\Phi\\left(\\frac{x - a}{\\sigma}\\right) - \\Phi\\left(\\frac{x - b}{\\sigma}\\right) \\right]
```

```julia
Bhattacharjee()        # equivalent to Bhattacharjee(0, 1, 1)
Bhattacharjee(σ)       # equivalent to Bhattacharjee(0, 1, σ)

params(d)        # Get the parameters, i.e. (a , b, σ)
```

External link:

*[Bhattacharjee Distributions on Technometrics](https://www.tandfonline.com/doi/abs/10.1080/00401706.1963.10490111)
"""
struct Bhattacharjee{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    a::T
    b::T
    σ::T
    Bhattacharjee{T}(a, b, σ) where {T<:Real} = new{T}(a, b, σ)
end

Bhattacharjee(a::Integer, b::Integer, σ::Integer; check_args::Bool=true) = Bhattacharjee(float(a), float(b), float(σ); check_args=check_args)

function Bhattacharjee(a::Real, b::Real, σ::Real; check_args::Bool=true)
    a, b, σ = promote(a, b, σ)
    @check_args Bhattacharjee (b, a < b) (σ, σ > zero(σ))
    return Bhattacharjee{typeof(a)}(a, b, σ)
end

Bhattacharjee() = Bhattacharjee{Float64}(0, 1, 1)
Bhattacharjee(σ) = Bhattacharjee{Float64}(0, 1, σ)

@distr_support Bhattacharjee -Inf Inf


# parameters

params(d::Bhattacharjee) = (d.a, d.b, d.σ)
@inline partype(d::Bhattacharjee{T}) where {T<:Real} = T

location(d::Bhattacharjee) = (d.a + d.b) / 2
scale(d::Bhattacharjee) = d.σ
shape(d::Bhattacharjee) = d.b - d.a

Base.eltype(::Type{Bhattacharjee{T}}) where {T} = T

#Statistic


Statistics.mean(d::Bhattacharjee) = (d.a + d.b)/2

Statistics.var(d::Bhattacharjee) = d.σ^2 + (d.b - d.a)^2 / 12

#### evaluate functions CDF, PDF, logPDF an CF

function Distributions.cdf(d::Bhattacharjee, x::Real)
    a, b, σ = d.a, d.b, d.σ
    t1, t2 = (x - a) / σ, (x - b) / σ
    (σ / (b - a)) * (t1 * Distributions.cdf(normal_dist, t1) - t2 * Distributions.cdf(normal_dist, t2) + Distributions.pdf(normal_dist, t1) - Distributions.pdf(normal_dist, t2))
end

function Distributions.pdf(d::Bhattacharjee, x::Real)
    a, b, σ = d.a, d.b, d.σ
    t1, t2 = (x - a) / σ, (x - b) / σ
    (1/(b - a)) * (Distributions.cdf(normal_dist, t1) - Distributions.cdf(normal_dist, t2))
end

function Distributions.logpdf(d::Bhattacharjee, x::Real)
    a, b, σ = d.a, d.b, d.σ
    t1, t2 = (x - a) / σ, (x - b) / σ 
    log(1/(b - a)) + log(Distributions.cdf(normal_dist, t1) - Distributions.cdf(normal_dist, t2))
end

function Distributions.quantile(d::Bhattacharjee, p::Real)
    a, b = d.a, d.b
    cdf_func(x) = Distributions.cdf(d, x) - p
    pdf_func(x) = Distributions.pdf(d, x)
    initial_point = a + p * (b - a) 
    q = Roots.find_zero((cdf_func, pdf_func), initial_point, Roots.Newton())
    return q 
end
function Distributions.cf(d::Bhattacharjee, t)
    a, b, σ = d.a, d.b, d.σ
    e_term = exp(-0.5 * σ^2 * t^2)
    u_term = (exp(1im * t * b) - exp(1im * t * a)) / (1im * t * (b - a))
    return e_term * u_term
end

function Distributions.mgf(d::Bhattacharjee, t::Real)
    a, b, σ = d.a, d.b, d.σ
    e_term = exp(0.5 * σ^2 * t^2)
    u_term = (exp(t * b) - exp(t * a)) / (t * (b - a))
    return e_term * u_term
end

## sampling 
function Distributions.rand(rng::Distributions.AbstractRNG, d::Bhattacharjee)
    a, b, σ = d.a, d.b, d.σ
    θ = rand(rng, Distributions.Uniform(a,b))
    return rand(rng, Distributions.Normal(θ,σ))
end