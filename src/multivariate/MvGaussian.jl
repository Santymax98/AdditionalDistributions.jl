"""
    MvGaussian(μ::AbstractVector, Σ::AbstractMatrix)

A multivariate Gaussian distribution backed by `Distributions.AbstractMvNormal`,
with rectangular cumulative probabilities evaluated by the custom
Genz/Richtmyer randomized quasi-Monte Carlo backend in
`AdditionalDistributions.jl`.

`MvGaussian` participates in the `Distributions.jl` `AbstractMvNormal`
interface, so standard functionality such as `pdf`, `logpdf`, `rand`,
`insupport`, `mode`, and `entropy` is inherited from `Distributions.jl`
through a small set of delegated primitives.

```julia
MvGaussian(Σ)        # zero-mean version
MvGaussian(μ, Σ)     # with explicit mean and covariance

params(d)            # (μ, Σ)
cdf(d, x)            # P(X₁ ≤ x₁, ..., Xₚ ≤ xₚ)
cdf(d, a, b)         # P(a ≤ X ≤ b)
cdf_result(d, a, b)  # structured numerical result
```
"""
struct MvGaussian{D<:Distributions.AbstractMvNormal} <: Distributions.AbstractMvNormal
    dist::D
end

MvGaussian(μ::AbstractVector, Σ::AbstractMatrix) =
    MvGaussian(Distributions.MvNormal(μ, Σ))

MvGaussian(Σ::AbstractMatrix) =
    MvGaussian(Distributions.MvNormal(Σ))


# ------------------------------------------------------------------
# AbstractMvNormal primitives
# ------------------------------------------------------------------

Base.length(d::MvGaussian) = length(d.dist)

Base.eltype(::Type{MvGaussian{D}}) where {D<:Distributions.AbstractMvNormal} =
    eltype(D)

Statistics.mean(d::MvGaussian) = Statistics.mean(d.dist)
Statistics.var(d::MvGaussian) = Statistics.var(d.dist)
Statistics.cov(d::MvGaussian) = Statistics.cov(d.dist)

params(d::MvGaussian) = params(d.dist)
partype(d::MvGaussian) = partype(d.dist)

Distributions.invcov(d::MvGaussian) =
    Distributions.invcov(d.dist)

Distributions.logdetcov(d::MvGaussian) =
    Distributions.logdetcov(d.dist)

Distributions.sqmahal(d::MvGaussian, x::AbstractVector) =
    Distributions.sqmahal(d.dist, x)

Distributions.sqmahal!(r::AbstractVector, d::MvGaussian, x::AbstractMatrix) =
    Distributions.sqmahal!(r, d.dist, x)

Distributions.gradlogpdf(d::MvGaussian, x::AbstractVector{<:Real}) =
    Distributions.gradlogpdf(d.dist, x)

Distributions._rand!(
    rng::Distributions.AbstractRNG,
    d::MvGaussian,
    x::AbstractVector,
) = Distributions._rand!(rng, d.dist, x)

Distributions._rand!(
    rng::Distributions.AbstractRNG,
    d::MvGaussian,
    x::AbstractMatrix,
) = Distributions._rand!(rng, d.dist, x)


# ------------------------------------------------------------------
# Rectangular Gaussian probabilities
# ------------------------------------------------------------------

"""
    cdf_result(d::Distributions.AbstractMvNormal, a, b;
               m=1000*length(a), abseps=1e-6, releps=1e-6,
               pivot=true, rng=Random.default_rng(),
               antithetic=false, batchsize=0, nshifts=10)

Estimate the rectangular probability `P(a ≤ X ≤ b)` for any
`Distributions.AbstractMvNormal` and return a [`CDFResult`](@ref).

This method is owned by `AdditionalDistributions.jl`, so it can legally extend
the CDF functionality of native `Distributions.MvNormal` objects without
extending `Distributions.cdf(::MvNormal, ...)`.
"""
function cdf_result(
    d::Distributions.AbstractMvNormal,
    a::AbstractVector,
    b::AbstractVector;
    m::Int=1000 * length(a),
    abseps::Real=1e-6,
    releps::Real=1e-6,
    pivot::Bool=true,
    rng=Random.default_rng(),
    antithetic::Bool=false,
    batchsize::Int=0,
    nshifts::Int=10,
)
    n = length(d)

    length(a) == n ||
        throw(DimensionMismatch("length(a) ≠ length(d)"))
    length(b) == n ||
        throw(DimensionMismatch("length(b) ≠ length(d)"))

    μ = Statistics.mean(d)
    Σ = Statistics.cov(d)

    size(Σ) == (n, n) ||
        throw(DimensionMismatch("Σ must be n×n"))

    T = promote_type(
        Float64,
        eltype(a),
        eltype(b),
        eltype(μ),
        eltype(Σ),
    )

    @inbounds for i in 1:n
        if b[i] < a[i]
            return CDFResult(
                zero(T),
                zero(T),
                0,
                0,
                :empty_rectangle,
            )
        end
    end

    if n == 1
        μ1 = T(μ[1])
        σ1 = sqrt(T(Σ[1, 1]))

        val =
            Distributions.cdf(
                Distributions.Normal(μ1, σ1),
                T(b[1]),
            ) -
            Distributions.cdf(
                Distributions.Normal(μ1, σ1),
                T(a[1]),
            )

        return CDFResult(
            clamp(T(val), zero(T), one(T)),
            zero(T),
            0,
            0,
            :univariate_exact,
        )
    end

    Σb = Matrix{T}(Σ)
    ab = T.(a)
    bb = T.(b)
    δb = T.(μ)

    res = mvtcdf(
        Σb,
        ab,
        bb;
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
        nshifts=nshifts,
    )

    return CDFResult(
        T(res.value),
        T(res.error),
        res.inform,
        m,
        :mvsort_rqmc,
    )
end


# Private scalar-value convenience API owned by AdditionalDistributions.
_cdf(
    d::Distributions.AbstractMvNormal,
    a::AbstractVector,
    b::AbstractVector;
    kwargs...,
) = cdf_result(d, a, b; kwargs...).value

function _cdf(
    d::Distributions.AbstractMvNormal,
    x::AbstractVector;
    kwargs...,
)
    length(x) == length(d) ||
        throw(DimensionMismatch("length(x) ≠ length(d)"))

    a = fill(-Inf, length(d))
    return _cdf(d, a, x; kwargs...)
end


# Public Distributions.cdf methods are defined only for our type.
function Distributions.cdf(
    d::MvGaussian,
    a::AbstractVector,
    b::AbstractVector;
    full::Bool=false,
    kwargs...,
)
    res = cdf_result(d, a, b; kwargs...)
    return full ? (res.value, res.error, res.inform) : res.value
end

function Distributions.cdf(
    d::MvGaussian,
    x::AbstractVector;
    full::Bool=false,
    kwargs...,
)
    length(x) == length(d) ||
        throw(DimensionMismatch("length(x) ≠ length(d)"))

    a = fill(-Inf, length(d))
    return Distributions.cdf(d, a, x; full=full, kwargs...)
end
