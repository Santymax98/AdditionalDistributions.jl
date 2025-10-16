struct MvTStudent{D<:Distributions.MvTDist} <: Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvTStudent(ν::Real, μ::AbstractVector, Σ::AbstractMatrix) = MvTStudent(MvTDist(ν, μ, Σ))

# ───────────────────────────────
#  Delegación de métodos estándar
# ───────────────────────────────
Statistics.mean(d::MvTStudent) = Statistics.mean(d.dist)
Statistics.cov(d::MvTStudent)  = Statistics.cov(d.dist)
Distributions.pdf(d::MvTStudent, x::AbstractVector) = Distributions.pdf(d.dist, x)
Distributions.logpdf(d::MvTStudent, x::AbstractVector) = Distributions.logpdf(d.dist, x)
Distributions.rand(rng::Distributions.AbstractRNG, d::MvTStudent) = Distributions.rand(rng, d.dist)
Distributions.insupport(d::MvTStudent, x::AbstractVector) = Distributions.insupport(d.dist, x)

# ───────────────────────────────
#  CDF por integración jerárquica
# ───────────────────────────────
function Distributions.cdf(d::MvTStudent, x::AbstractVector;
        r::Int=12, m::Int=10_000, rng=Random.default_rng(),
        batch::Int=96, cache_kib::Int=256, full::Bool=false)

    # Extrae parámetros desde tu wrapper (campos según tu MvTDist)
    μ  = d.dist.μ
    Σ  = d.dist.Σ
    νr = d.dist.df               # en tu tipo: df (no ν)

    # Centra y llama al núcleo t
    Δ = x .- μ
    TΣ = eltype(Σ)
    ν  = TΣ(νr)                  # asegura el tipo numérico consistente

    val, err = _cdf_t_base(Σ, Δ; ν=ν, r=r, m=m, rng=rng, batch=batch, cache_kib=cache_kib)
    return full ? (val, err) : val
end
