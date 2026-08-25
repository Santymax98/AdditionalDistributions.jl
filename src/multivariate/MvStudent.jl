"""
    MvTStudent(ν::Real, μ::AbstractVector, Σ::AbstractMatrix)

A multivariate Student's t distribution backed by a native
`Distributions.AbstractMvTDist`, with rectangular cumulative probabilities
evaluated by the custom Genz-Bretz/Richtmyer randomized QMC backend in
`AdditionalDistributions.jl`.

`Σ` is the scale/scatter matrix used by `Distributions.MvTDist`; it is not
generally the covariance matrix.

```julia
MvTStudent(ν, Σ)
MvTStudent(ν, μ, Σ)

params(d)            # (ν, μ, Σ)
cdf(d, x)            # P(X₁ ≤ x₁, ..., Xₚ ≤ xₚ)
cdf(d, a, b)         # P(a ≤ X ≤ b)
cdf_result(d, a, b)  # structured numerical result
```
"""
struct MvTStudent{D<:Distributions.AbstractMvTDist} <:
       Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvTStudent(
    ν::Real,
    μ::AbstractVector,
    Σ::AbstractMatrix,
) = MvTStudent(Distributions.MvTDist(ν, collect(μ), Matrix(Σ)))

MvTStudent(
    ν::Real,
    Σ::AbstractMatrix,
) = MvTStudent(
    Distributions.MvTDist(
        ν,
        zeros(eltype(Σ), size(Σ, 1)),
        Matrix(Σ),
    ),
)


# ------------------------------------------------------------------
# Standard wrapper delegation
# ------------------------------------------------------------------

Base.length(d::MvTStudent) = length(d.dist)

Base.eltype(::Type{MvTStudent{D}}) where {D<:Distributions.AbstractMvTDist} =
    eltype(D)

Statistics.mean(d::MvTStudent) = Statistics.mean(d.dist)
Statistics.var(d::MvTStudent) = Statistics.var(d.dist)
Statistics.cov(d::MvTStudent) = Statistics.cov(d.dist)

Distributions.pdf(d::MvTStudent, x::AbstractVector) =
    Distributions.pdf(d.dist, x)

Distributions.logpdf(d::MvTStudent, x::AbstractVector) =
    Distributions.logpdf(d.dist, x)

Distributions.rand(
    rng::Distributions.AbstractRNG,
    d::MvTStudent,
) = Distributions.rand(rng, d.dist)

Distributions.rand(
    rng::Distributions.AbstractRNG,
    d::MvTStudent,
    n::Int,
) = Distributions.rand(rng, d.dist, n)

Distributions.insupport(d::MvTStudent, x::AbstractVector) =
    Distributions.insupport(d.dist, x)

Distributions.scale(d::MvTStudent) =
    Distributions.scale(d.dist)

Distributions.invcov(d::MvTStudent) =
    Distributions.invcov(d.dist)

StatsBase.mode(d::MvTStudent) =
    StatsBase.mode(d.dist)

StatsBase.entropy(d::MvTStudent) =
    StatsBase.entropy(d.dist)

params(d::MvTStudent) = params(d.dist)
partype(d::MvTStudent) = partype(d.dist)


# ------------------------------------------------------------------
# Rectangular Student-t probabilities
# ------------------------------------------------------------------

"""
    cdf_result(d::Distributions.AbstractMvTDist, a, b;
               m=max(100_000, 10_000*length(a)),
               abseps=1e-6, releps=1e-6, pivot=true,
               rng=Random.default_rng(), antithetic=false,
               batchsize=0, nshifts=8)

Estimate `P(a ≤ X ≤ b)` for a native `Distributions.AbstractMvTDist` and
return a [`CDFResult`](@ref).

The calculation uses the t distribution's location `μ` and scale/scatter
matrix `Σ`, not its statistical mean and covariance.
"""
function cdf_result(
    d::Distributions.AbstractMvTDist,
    a::AbstractVector,
    b::AbstractVector;
    m::Int=max(100_000, 10_000 * length(a)),
    abseps::Real=1e-6,
    releps::Real=1e-6,
    pivot::Bool=true,
    rng=Random.default_rng(),
    antithetic::Bool=false,
    batchsize::Int=0,
    nshifts::Int=8,
)
    n = length(d)

    length(a) == n ||
        throw(DimensionMismatch("length(a) ≠ length(d)"))
    length(b) == n ||
        throw(DimensionMismatch("length(b) ≠ length(d)"))

    # Distributions.jl's AbstractMvTDist interface itself relies on these
    # canonical fields. Importantly, μ is location and Σ is scale/scatter.
    μ = d.μ
    Σscale = d.Σ
    ν = d.df

    size(Σscale) == (n, n) ||
        throw(DimensionMismatch("Σ must be n×n"))

    T = promote_type(
        Float64,
        eltype(a),
        eltype(b),
        eltype(μ),
        eltype(Σscale),
        typeof(ν),
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
        σ1 = sqrt(T(Σscale[1, 1]))
        td = Distributions.TDist(T(ν))

        za = (T(a[1]) - μ1) / σ1
        zb = (T(b[1]) - μ1) / σ1

        val = Distributions.cdf(td, zb) - Distributions.cdf(td, za)

        return CDFResult(
            clamp(T(val), zero(T), one(T)),
            zero(T),
            0,
            0,
            :univariate_exact_t,
        )
    end

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

            for i in (j + 1):n
                cij =
                    T(Σscale[i, j]) *
                    invD[i] *
                    invD[j]

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

    res = mvtcdf(
        C,
        aS,
        bS;
        ν=T(ν),
        δ=δS,
        maxpts=m,
        abseps=abseps,
        releps=releps,
        assume_correlation=true,
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
        :mvsort_rqmc_t,
    )
end


cdf_result(
    d::MvTStudent,
    a::AbstractVector,
    b::AbstractVector;
    kwargs...,
) = cdf_result(d.dist, a, b; kwargs...)


# Private scalar-value convenience API owned by AdditionalDistributions.
_cdf(
    d::Distributions.AbstractMvTDist,
    a::AbstractVector,
    b::AbstractVector;
    kwargs...,
) = cdf_result(d, a, b; kwargs...).value

function _cdf(
    d::Distributions.AbstractMvTDist,
    x::AbstractVector;
    kwargs...,
)
    length(x) == length(d) ||
        throw(DimensionMismatch("length(x) ≠ length(d)"))

    a = fill(-Inf, length(d))
    return _cdf(d, a, x; kwargs...)
end


# Public Distributions.cdf methods are defined only for our wrapper type.
function Distributions.cdf(
    d::MvTStudent,
    a::AbstractVector,
    b::AbstractVector;
    full::Bool=false,
    kwargs...,
)
    res = cdf_result(d, a, b; kwargs...)
    return full ? (res.value, res.error, res.inform) : res.value
end

function Distributions.cdf(
    d::MvTStudent,
    x::AbstractVector;
    full::Bool=false,
    kwargs...,
)
    length(x) == length(d) ||
        throw(DimensionMismatch("length(x) ≠ length(d)"))

    a = fill(-Inf, length(d))
    return Distributions.cdf(d, a, x; full=full, kwargs...)
end
