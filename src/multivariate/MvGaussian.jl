"""
    MvGaussian(μ::AbstractVector, Σ::AbstractMatrix)

A *Multivariate Gaussian (Normal) distribution* equivalent in behavior to
[`Distributions.MvNormal`](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvNormal),
but using a **custom folded batch randomized quasi-Monte Carlo integrator** for the cumulative
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
Distributions.rand(rng::Distributions.AbstractRNG, d::MvGaussian, n::Int) = Distributions.rand(rng, d.dist, n)
Distributions.insupport(d::MvGaussian, x::AbstractVector)         = Distributions.insupport(d.dist, x)
params(d::MvGaussian) = (Statistics.mean(d), Statistics.cov(d))
Base.length(d::MvGaussian) = length(d.dist)


"""
    cdf_result(d::MvGaussian, a, b; m=1000*length(a), abseps=1e-6,
               releps=1e-6, pivot=true, rng=Random.default_rng(),
               antithetic=false, batchsize=0, nshifts=12)

Estimate the rectangular probability `P(a ≤ X ≤ b)` for a multivariate
Gaussian distribution and return a [`CDFResult`](@ref).

The Gaussian path uses MVSORT reordering, a Genz-style conditional
transformation, folded randomized Richtmyer QMC points, and batch evaluation.
If the covariance matrix is diagonal, the probability is evaluated exactly by
factorization.

Use `cdf(d, a, b)` for the scalar probability, or `cdf(d, a, b; full=true)` for
the legacy `(value, error, inform)` tuple.
"""
function cdf_result(d::MvGaussian,
                    a::AbstractVector, b::AbstractVector;
                    m::Int=1000*length(a),
                    abseps::Real=1e-6,
                    releps::Real=1e-6,
                    pivot::Bool=true,
                    rng=Random.default_rng(),
                    antithetic::Bool=false,
                    batchsize::Int=0,
                    nshifts::Int=12)

    n = length(a)
    length(b) == n || throw(DimensionMismatch("length(a) ≠ length(b)"))

    μ = Statistics.mean(d)
    Σ = Statistics.cov(d)

    size(Σ) == (n, n) || throw(DimensionMismatch("Σ must be n×n"))

    T = promote_type(Float64, eltype(a), eltype(b), eltype(μ), eltype(Σ))

    if any(b .< a)
        return CDFResult(zero(T), zero(T), 0, 0, :empty_rectangle)
    end

    if n == 1
        μ1 = μ[1]
        σ1 = sqrt(Σ[1, 1])
        val = Distributions.cdf(Distributions.Normal(μ1, σ1), b[1]) -
              Distributions.cdf(Distributions.Normal(μ1, σ1), a[1])
        val = clamp(T(val), zero(T), one(T))
        return CDFResult(val, zero(T), 0, 0, :univariate_exact)
    end

    Σb = T.(Σ)
    ab = T.(a)
    bb = T.(b)
    δb = T.(μ)

    res = mvtcdf(Σb, ab, bb;
                 ν=0,
                 δ=δb,
                 maxpts=m,
                 abseps=abseps,
                 releps=releps,
                 assume_correlation=false,
                 pivot=pivot,
                 antithetic=antithetic,
                 rng=rng,
                 batchsize=batchsize,
                 nshifts=nshifts)

    return CDFResult(T(res.value), T(res.error), res.inform, m, :mvsort_rqmc)
end

function Distributions.cdf(d::MvGaussian,
                           a::AbstractVector, b::AbstractVector;
                           full::Bool=false,
                           kwargs...)

    res = cdf_result(d, a, b; kwargs...)
    return full ? (res.value, res.error, res.inform) : res.value
end