# ─────────────────────────────────────────────────────────────
# 1. Tipo envoltorio (Wrapper sobre MvNormal)
# ─────────────────────────────────────────────────────────────
struct MvGaussian{D<:Distributions.MvNormal} <: Distributions.ContinuousMultivariateDistribution
    dist::D
end

MvGaussian(μ::AbstractVector, Σ::AbstractMatrix) = MvGaussian(Distributions.MvNormal(μ, Σ))

# ─────────────────────────────────────────────────────────────
# 2. Delegaciones
# ─────────────────────────────────────────────────────────────
Statistics.mean(d::MvGaussian) = Statistics.mean(d.dist)
Statistics.cov(d::MvGaussian)  = Statistics.cov(d.dist)
Distributions.pdf(d::MvGaussian, x::AbstractVector) = Distributions.pdf(d.dist, x)
Distributions.logpdf(d::MvGaussian, x::AbstractVector) = Distributions.logpdf(d.dist, x)
Distributions.rand(rng::Distributions.AbstractRNG, d::MvGaussian) = Distributions.rand(rng, d.dist)

Distributions.insupport(d::MvGaussian, x::AbstractVector) = Distributions.insupport(d.dist, x)

# ─────────────────────────────────────────────────────────────
# 3. Utilidades gaussianas
# ─────────────────────────────────────────────────────────────
@inline Φ(x::T) where {T<:AbstractFloat} = T(0.5) * SpecialFunctions.erfc(-x / sqrt(T(2)))
@inline Φinv(p::T) where {T<:AbstractFloat} = sqrt(T(2)) * SpecialFunctions.erfcinv(T(2) * (T(1) - p))

# ─────────────────────────────────────────────────────────────
# 4. CDF QMC (tu implementación optimizada)
# ─────────────────────────────────────────────────────────────
function _cdf(Σ::AbstractMatrix{T}, x::AbstractVector{T};
              r::Int=12, m::Int=10_000, rng=Random.default_rng(),
              batch::Int=96, cache_kib::Int=256) where {T<:AbstractFloat}
    d = length(x)
    nv = max(m ÷ r, 1)
    B = if batch ≤ 0
        bytes_row = 2d * sizeof(T)
        B_auto = max(64, cld(cache_kib * 1024, bytes_row))
        min(nv, B_auto)
    else
        min(nv, batch)
    end
    σ = sqrt.(LinearAlgebra.diag(Σ))
    R = Σ ./ (σ * σ')
    b0 = x ./ σ
    P = sortperm(Φ.(b0))
    R1 = @view R[P, P]

    F = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(R1), LinearAlgebra.RowMaximum(); check=false)
    L = Matrix(F.L)
    b = b0[P][F.p]

    yB, uB = Matrix{T}(undef, B, d), Matrix{T}(undef, B, d)
    tv, βv, pv = (Vector{T}(undef, B) for _ in 1:3)
    tmpd, shift = Vector{T}(undef, d), Vector{T}(undef, d)
    sob = Sobol.SobolSeq(d)

    ε, zero_T, one_T = eps(T), zero(T), one(T)
    tiny = floatmin(T)
    β1 = Φ(b[1] / L[1,1])

    p_mean, m2 = zero_T, zero_T

    @inbounds for j in 1:r
        Distributions.rand!(rng, shift)
        sum_pv = zero_T

        for bstart in 1:B:nv
            cnt = min(B, nv - bstart + 1)

            for i in 1:cnt
                Sobol.next!(sob, tmpd)
                @simd for k in 1:d
                    u = tmpd[k] + shift[k]
                    uB[i,k] = u - (u ≥ one_T ? one_T : zero_T)
                end
            end

            if β1 > ε
                @simd for i in 1:cnt
                    yB[i,1] = Φinv(clamp(uB[i,1] * β1, tiny, one_T - ε))
                    pv[i] = β1
                end
            else
                @simd for i in 1:cnt
                    pv[i] = zero_T
                end
            end

            for k in 2:d
                Lkk_inv = inv(L[k,k])

                for i in 1:cnt
                    s = zero_T
                    @simd for j2 in 1:(k-1)
                        s += L[k,j2] * yB[i,j2]
                    end
                    tv[i] = (b[k] - s) * Lkk_inv
                end

                @simd for i in 1:cnt
                    βv[i] = Φ(tv[i])
                end

                for i in 1:cnt
                    β = βv[i]
                    if (pv[i] > zero_T) & (β > ε)
                        yB[i,k] = Φinv(clamp(uB[i,k] * β, tiny, one_T - ε))
                        pv[i] *= β
                    else
                        pv[i] = zero_T
                    end
                end
            end

            sblk = zero_T
            @simd for i in 1:cnt
                sblk += pv[i]
            end
            sum_pv += sblk
        end

        mean_j = sum_pv / nv
        old_mean = p_mean
        p_mean += (mean_j - p_mean) / j
        m2 += (mean_j - old_mean) * (mean_j - p_mean)
    end

    # ERROR ESTIMADO COMPARABLE (3 errores estándar)
    err = r > 1 ? sqrt((m2 / (r - 1)) / r) : T(NaN)
    
    return (p_mean, err)
end

# ─────────────────────────────────────────────────────────────
# 5. CDF despachado para MvGaussian
# ─────────────────────────────────────────────────────────────
function Distributions.cdf( d::MvGaussian, x::AbstractVector; r::Int=12, m::Int=10_000, 
                rng=Random.default_rng(), batch::Int=96, cache_kib::Int=256, full::Bool=false)
    μ = Statistics.mean(d.dist)
    Σ = Statistics.cov(d.dist)
    Δ = x .- μ
    val, err = _cdf(Σ, Δ; r=r, m=m, rng=rng, batch=batch, cache_kib=cache_kib)
    return full ? (val, err) : val
end
