struct MvGaussian{D<:Distributions.MvNormal} <: Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvGaussian(μ::AbstractVector, Σ::AbstractMatrix) = MvGaussian(Distributions.MvNormal(μ, Σ))

# ───────────────────────────────
#  Delegación de métodos estándar
# ───────────────────────────────
Statistics.mean(d::MvGaussian) = Statistics.mean(d.dist)
Statistics.cov(d::MvGaussian)  = Statistics.cov(d.dist)
Distributions.pdf(d::MvGaussian, x::AbstractVector) = Distributions.pdf(d.dist, x)
Distributions.logpdf(d::MvGaussian, x::AbstractVector) = Distributions.logpdf(d.dist, x)
Distributions.rand(rng::Distributions.AbstractRNG, d::MvGaussian) = Distributions.rand(rng, d.dist)

Distributions.insupport(d::MvGaussian, x::AbstractVector) = Distributions.insupport(d.dist, x)

# ─────────────────────────────────────────────────────────────
# 5. CDF despachado para MvGaussian
# ─────────────────────────────────────────────────────────────

function Distributions.cdf(d::MvGaussian, x::AbstractVector; r::Int=12, m::Int=10_000, 
        rng=Random.default_rng(), batch::Int=96, cache_kib::Int=256, full::Bool=false)
    μ = Statistics.mean(d.dist)
    Σ = Statistics.cov(d.dist)
    Δ = x .- μ
    val, err = _cdf_base(Σ, Δ; r=r, m=m, rng=rng, batch=batch, cache_kib=cache_kib)
    return full ? (val, err) : val
end