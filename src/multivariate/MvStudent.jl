"""
    MvTStudent(ν::Real, μ::AbstractVector, Σ::AbstractMatrix)

A *Multivariate Student's t distribution* equivalent in behavior to
[`Distributions.MvTDist`](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvTDist),
but using a **custom Quasi-Monte Carlo (QMC) integrator** for the
rectangular cumulative distribution function (`cdf`).

This type preserves all the standard functionality of `MvTDist` — including
`pdf`, `logpdf`, `rand`, `mean`, and `cov` — while providing its own
implementation of `cdf(a, b)` for rectangular probabilities under a
multivariate t law.

```julia
MvTStudent(ν, Σ)             # zero-mean version with df=ν
MvTStudent(ν, μ, Σ)          # with degrees of freedom ν, mean μ, and covariance Σ

params(d)                    # returns (ν, μ, Σ)
cdf(d, a, b)                 # evaluates P(a ≤ X ≤ b)
```

External link:

* [Multivariate Student's t-distribution - Wikipedia](https://en.wikipedia.org/wiki/Multivariate_t-distribution)
"""
struct MvTStudent{D<:Distributions.MvTDist} <: Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvTStudent(ν::Real, μ::AbstractVector, Σ::AbstractMatrix) = MvTStudent(Distributions.MvTDist(ν, μ, Σ))
MvTStudent(ν::Real, Σ::AbstractMatrix) = MvTStudent(Distributions.MvTDist(ν, zeros(size(Σ, 1)), Σ))


# ──────────────────────────
# Standard Method Delegation
# ──────────────────────────
Statistics.mean(d::MvTStudent) = Statistics.mean(d.dist)
Statistics.cov(d::MvTStudent)  = Statistics.cov(d.dist)
Distributions.pdf(d::MvTStudent, x::AbstractVector) = Distributions.pdf(d.dist, x)
Distributions.logpdf(d::MvTStudent, x::AbstractVector) = Distributions.logpdf(d.dist, x)
Distributions.rand(rng::Distributions.AbstractRNG, d::MvTStudent) = Distributions.rand(rng, d.dist)
Distributions.insupport(d::MvTStudent, x::AbstractVector) = Distributions.insupport(d.dist, x)
params(d::MvTStudent) = (d.dist.df, d.dist.μ, d.dist.Σ)

# ───────────────────────────────
# CDF by GENZ-BRETZ
# ───────────────────────────────
function Distributions.cdf(d::MvTStudent,
                           a::AbstractVector, b::AbstractVector;
                           m::Int = 1000*length(a),
                           abseps::Real = 1e-6,
                           releps::Real = 1e-6,
                           pivot::Bool = true,
                           rng = Random.default_rng(),
                           full::Bool = false,
                           antithetic::Bool=false)

    n = length(a); @assert length(b) == n
    μ = d.dist.μ                       # location, also defined when ν ≤ 1
    Σscale = d.dist.Σ                   # scatter matrix, not covariance
    ν = d.dist.df

    @assert size(Σscale) == (n,n)
    @assert all(b .>= a)

    T = promote_type(eltype(Σscale), eltype(a), eltype(b), eltype(μ))
    Ddiag = sqrt.(LinearAlgebra.diag(Σscale))
    invD  = 1 ./(T.(Ddiag))
    C     = LinearAlgebra.Symmetric(LinearAlgebra.Diagonal(invD) * T.(Σscale) * LinearAlgebra.Diagonal(invD))
    aS    = (T.(a) .- T.(μ)) .* invD
    bS    = (T.(b) .- T.(μ)) .* invD
    δS    = zeros(T, n)

    res = mvtcdf(C, aS, bS; ν=ν, δ=δS,
                 maxpts=m, abseps=abseps, releps=releps,
                 assume_correlation=true, pivot=pivot,
                 antithetic=antithetic, rng=rng)

    return full ? (res.value, res.error, res.inform) : res.value
end


