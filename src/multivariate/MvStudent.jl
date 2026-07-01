"""
    MvTStudent(ν::Real, μ::AbstractVector, Σ::AbstractMatrix)

A *Multivariate Student's t distribution* equivalent in behavior to
[`Distributions.MvTDist`](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvTDist),
but using a custom pure-Julia Genz-Bretz/Richtmyer randomized QMC integrator
for rectangular cumulative distribution functions.

This type preserves the standard functionality of `MvTDist` — including
`pdf`, `logpdf`, `rand`, `mean`, and `cov` — while providing its own
implementation of `cdf(a, b)` for rectangular probabilities under a
multivariate t law.

```julia
MvTStudent(ν, Σ)             # zero-location version with df=ν
MvTStudent(ν, μ, Σ)          # with degrees of freedom ν, location μ, and scale Σ

params(d)                    # returns (ν, μ, Σ)
cdf(d, a, b)                 # evaluates P(a ≤ X ≤ b)
cdf_result(d, a, b)          # returns value, error, inform and algorithm metadata
```
"""
struct MvTStudent{D<:Distributions.MvTDist} <: Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvTStudent(ν::Real, μ::AbstractVector, Σ::AbstractMatrix) = MvTStudent(Distributions.MvTDist(ν, μ, Σ))
MvTStudent(ν::Real, Σ::AbstractMatrix) = MvTStudent(Distributions.MvTDist(ν, zeros(size(Σ, 1)), Σ))

# ──────────────────────────
# Standard method delegation
# ──────────────────────────
Statistics.mean(d::MvTStudent) = Statistics.mean(d.dist)
Statistics.cov(d::MvTStudent) = Statistics.cov(d.dist)
Distributions.pdf(d::MvTStudent, x::AbstractVector) = Distributions.pdf(d.dist, x)
Distributions.logpdf(d::MvTStudent, x::AbstractVector) = Distributions.logpdf(d.dist, x)
Distributions.rand(rng::Distributions.AbstractRNG, d::MvTStudent) = Distributions.rand(rng, d.dist)
Distributions.insupport(d::MvTStudent, x::AbstractVector) = Distributions.insupport(d.dist, x)
params(d::MvTStudent) = (d.dist.df, d.dist.μ, d.dist.Σ)

# ───────────────────────────────
# Rectangular CDF by Genz-Bretz/Richtmyer RQMC
# ───────────────────────────────
function cdf_result(d::MvTStudent,
                    a::AbstractVector, b::AbstractVector;
                    m::Int=max(100_000, 10_000 * length(a)),
                    abseps::Real=1e-6,
                    releps::Real=1e-6,
                    pivot::Bool=true,
                    rng=Random.default_rng(),
                    antithetic::Bool=false,
                    batchsize::Int=0,
                    nshifts::Int=16)

    n = length(a)
    length(b) == n || throw(DimensionMismatch("length(a) ≠ length(b)"))

    μ = d.dist.μ
    Σscale = d.dist.Σ
    ν = d.dist.df

    size(Σscale) == (n, n) || throw(DimensionMismatch("Σ must be n×n"))

    T = promote_type(Float64, eltype(a), eltype(b), eltype(μ), eltype(Σscale))

    @inbounds for i in 1:n
        if b[i] < a[i]
            return CDFResult(zero(T), zero(T), 0, 0, :empty_rectangle)
        end
    end

    # Standardize to zero-location, unit diagonal scale. For the multivariate t,
    # Σscale is the scatter/scale matrix used by MvTDist, not necessarily the
    # covariance matrix. This diagonal scaling preserves the t dependence and
    # passes a correlation-like scale matrix to the core integrator.
    Ddiag = Vector{T}(undef, n)
    invD = Vector{T}(undef, n)
    @inbounds for i in 1:n
        Ddiag[i] = sqrt(T(Σscale[i, i]))
        invD[i] = inv(Ddiag[i])
    end

    C = Matrix{T}(undef, n, n)
    @inbounds begin
        for j in 1:n
            C[j, j] = one(T)
            for i in j+1:n
                cij = T(Σscale[i, j]) * invD[i] * invD[j]
                C[i, j] = cij
                C[j, i] = cij
            end
        end
    end

    aS = Vector{T}(undef, n)
    bS = Vector{T}(undef, n)
    δS = zeros(T, n)
    @inbounds for i in 1:n
        aS[i] = (T(a[i]) - T(μ[i])) * invD[i]
        bS[i] = (T(b[i]) - T(μ[i])) * invD[i]
    end

    res = mvtcdf(C, aS, bS;
                 ν=ν,
                 δ=δS,
                 maxpts=m,
                 abseps=abseps,
                 releps=releps,
                 assume_correlation=true,
                 pivot=pivot,
                 antithetic=antithetic,
                 rng=rng,
                 batchsize=batchsize,
                 nshifts=nshifts)

    return CDFResult(T(res.value), T(res.error), res.inform, m, :mvsort_rqmc_t)
end

function Distributions.cdf(d::MvTStudent,
                           a::AbstractVector, b::AbstractVector;
                           full::Bool=false,
                           kwargs...)

    res = cdf_result(d, a, b; kwargs...)
    return full ? (res.value, res.error, res.inform) : res.value
end
