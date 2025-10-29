"""
    MvGaussian(μ::AbstractVector, Σ::AbstractMatrix)

A *Multivariate Gaussian (Normal) distribution* equivalent in behavior to
[`Distributions.MvNormal`](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvNormal),
but using a **custom Quasi-Monte Carlo (QMC) integrator** for the cumulative
distribution function (`cdf`).

This type preserves all the standard functionality of `MvNormal` — including
`pdf`, `logpdf`, `rand`, `mean`, and `cov` — while providing its own
implementation of `cdf(a, b)` for rectangular probabilities under a
multivariate normal law.

```julia
MvGaussian(Σ)        # zero-mean version
MvGaussian(μ, Σ)     # with explicit mean and covariance
params(d)            # returns (μ, Σ)
cdf(d, a, b)         # evaluates P(a ≤ X ≤ b)
```

External link:

* [Multivariate normal distribution - Wikipedia](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
"""
struct MvGaussian{D<:Distributions.MvNormal} <: Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvGaussian(μ::AbstractVector, Σ::AbstractMatrix) = MvGaussian(Distributions.MvNormal(μ, Σ))
MvGaussian(Σ::AbstractMatrix) = MvGaussian(Distributions.MvNormal(zeros(size(Σ, 1)), Σ))


# Delegaciones
Statistics.mean(d::MvGaussian) = Statistics.mean(d.dist)
Statistics.cov(d::MvGaussian)  = Statistics.cov(d.dist)
Distributions.pdf(d::MvGaussian, x::AbstractVector)    = Distributions.pdf(d.dist, x)
Distributions.logpdf(d::MvGaussian, x::AbstractVector) = Distributions.logpdf(d.dist, x)
Distributions.rand(rng::Distributions.AbstractRNG, d::MvGaussian) = Distributions.rand(rng, d.dist)
Distributions.insupport(d::MvGaussian, x::AbstractVector)         = Distributions.insupport(d.dist, x)

# MvGaussian.jl
function Distributions.cdf(d::MvGaussian,
                           a::AbstractVector, b::AbstractVector;
                           m::Int=1000*length(a),
                           abseps::Real=1e-6, releps::Real=1e-6,
                           pivot::Bool=true,
                           rng=Random.default_rng(),
                           full::Bool=false,
                           antithetic::Bool=true)

    n = length(a)
    @assert length(b) == n "length(a) ≠ length(b)"
    Σ = Statistics.cov(d)
    @assert size(Σ) == (n,n) "Σ no es cuadrada de tamaño n×n"
    if n == 1
        throw(ErrorException("Σ dimension 1 not supported"))
    end
    if any(b .< a)
        throw(ArgumentError("se encontró a[i] > b[i]"))
    end

    μ = Statistics.mean(d)
    δ = μ

    Tb = promote_type(eltype(Σ), eltype(a), eltype(b), eltype(δ))
    Σb = Tb.(Σ); ab = Tb.(a); bb = Tb.(b); δb = Tb.(δ)

    res = mvtcdf(Σb, ab, bb; ν=0, δ=δb,
                 maxpts=m, abseps=abseps, releps=releps,
                 assume_correlation=false, pivot=pivot,
                 antithetic=antithetic, rng=rng)

    return full ? (res.value, res.error, res.inform) : res.value
end