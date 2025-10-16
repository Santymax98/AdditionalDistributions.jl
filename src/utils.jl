@inline Φ(x::T) where {T<:AbstractFloat} = @fastmath T(0.5) * SpecialFunctions.erfc(-x / sqrt(T(2)))
@inline Φinv(p::T) where {T<:AbstractFloat} = @fastmath sqrt(T(2)) * SpecialFunctions.erfcinv(T(2) * (T(1) - p))

function _cdf_base(Σ::AbstractMatrix{T}, x::AbstractVector{T};
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

    err = r > 1 ? sqrt((m2 / (r - 1)) / r) : T(NaN)
    return (p_mean, err)
end



function _cdf_t_base(Σ::AbstractMatrix{T}, x::AbstractVector{T};
                     ν::T, r::Int=12, m::Int=10_000, rng=Random.default_rng(),
                     batch::Int=96, cache_kib::Int=256) where {T<:AbstractFloat}

    if !isfinite(ν)
        return _cdf_base(Σ, x; r=r, m=m, rng=rng, batch=batch, cache_kib=cache_kib)
    end

    d  = length(x)
    nv = max(m ÷ r, 1)
    B  = if batch ≤ 0
        bytes_row = 2d * sizeof(T)
        B_auto = max(64, cld(cache_kib * 1024, bytes_row))
        min(nv, B_auto)
    else
        min(nv, batch)
    end

    σ  = sqrt.(LinearAlgebra.diag(Σ))
    R  = Σ ./ (σ * σ')
    b0 = x ./ σ

    P  = sortperm(Φ.(b0))  # w>0 no cambia el orden
    R1 = @view R[P, P]

    F  = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(R1), LinearAlgebra.RowMaximum(); check=false)
    L  = Matrix(F.L)
    b  = b0[P][F.p]

    yB  = Matrix{T}(undef, B, d)
    uB  = Matrix{T}(undef, B, d + 1)   # +1 para u_s
    tv  = Vector{T}(undef, B)
    βv  = Vector{T}(undef, B)
    pv  = Vector{T}(undef, B)
    wv  = Vector{T}(undef, B)
    tmp = Vector{T}(undef, d + 1)
    shift = Vector{T}(undef, d + 1)
    sob = Sobol.SobolSeq(d + 1)

    ε, zero_T, one_T = eps(T), zero(T), one(T)
    tiny = floatmin(T)
    p_mean, m2 = zero_T, zero_T

    a = T(0.5) * ν  # parámetro de la gamma incompleta

    @inbounds for j in 1:r
        Distributions.rand!(rng, shift)
        sum_pv = zero_T

        for bstart in 1:B:nv
            cnt = min(B, nv - bstart + 1)

            # Sobol(d+1) + corrimiento
            for i in 1:cnt
                Sobol.next!(sob, tmp)
                @simd for k in 1:(d+1)
                    u = tmp[k] + shift[k]
                    uB[i,k] = u - (u ≥ one_T ? one_T : zero_T)
                end
            end

            # Mezcla: s = 2 * gamma_inc_inv(a, p=u, q=1-u); w = √(s/ν)
            @simd for i in 1:cnt
                us   = clamp(uB[i,1], tiny, one_T - ε)
                xgam = SpecialFunctions.gamma_inc_inv(a, us, one_T - us)  # devuelve x tal que P(a,x)=us
                s    = 2*xgam
                wv[i] = sqrt(s / ν)
            end

            # Primer margen condicional (dependiente de w_i)
            L11 = L[1,1]
            @simd for i in 1:cnt
                β = Φ((wv[i]*b[1]) / L11)
                if β > ε
                    pv[i]   = β
                    u2      = clamp(uB[i,2] * β, tiny, one_T - ε)
                    yB[i,1] = Φinv(u2)
                else
                    pv[i]   = zero_T
                    yB[i,1] = zero_T
                end
            end

            # Márgenes k = 2..d
            for k in 2:d
                Lkk    = L[k,k]
                invLkk = one_T / Lkk
                for i in 1:cnt
                    ssum = zero_T
                    @simd for j2 in 1:(k-1)
                        ssum += L[k,j2] * yB[i,j2]
                    end
                    tv[i] = (wv[i]*b[k] - ssum) * invLkk
                end

                @simd for i in 1:cnt
                    βv[i] = Φ(tv[i])
                end

                for i in 1:cnt
                    β = βv[i]
                    if (pv[i] > zero_T) & (β > ε)
                        uk = clamp(uB[i, k+1] * β, tiny, one_T - ε)  # usa columna k+1
                        yB[i,k] = Φinv(uk)
                        pv[i]  *= β
                    else
                        pv[i] = zero_T
                        yB[i,k] = zero_T
                    end
                end
            end

            sblk = zero_T
            @simd for i in 1:cnt
                sblk += pv[i]
            end
            sum_pv += sblk
        end

        mean_j  = sum_pv / nv
        δ       = mean_j - p_mean
        p_mean += δ / j
        m2     += δ * (mean_j - p_mean)
    end

    err = r > 1 ? sqrt((m2 / (r - 1)) / r) : T(NaN)
    return (p_mean, err)
end
