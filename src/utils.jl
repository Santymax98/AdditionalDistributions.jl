# ─────────────────────────────────────────────────────────────
# Utils Normal(0,1)
# ─────────────────────────────────────────────────────────────
@inline _Z(T::Type{<:Real}) = Distributions.Normal(zero(T), one(T))
@inline _cdf(x::T) where {T<:Real} = Distributions.cdf(_Z(T), x)
@inline _pdf(x::T) where {T<:Real} = Distributions.pdf(_Z(T), x)
@inline _qf_std(p::T) where {T<:Real} = Distributions.quantile(_Z(T), p)

# χ² → scale for t: R = sqrt(χ²_ν)/sqrt(ν)
@inline function _scale_t(ν::Real, u::Real)
    @assert ν > 0
    return sqrt(Distributions.quantile(Distributions.Chisq(ν), u)) / sqrt(ν)
end

# ────────────────────────
# Lower triangular packing
# ────────────────────────
struct PackedMatrix{T<:Real}
    data::Vector{T}
    n::Int
end
# 1-based index for (i,j)
@inline _lin(i::Int, j::Int) = i ≥ j ? (i*(i-1))÷2 + j : (j*(j-1))÷2 + i
@propagate_inbounds @inline function _get_cov(C::PackedMatrix, i::Int, j::Int)
    return C.data[_lin(i,j)]
end

# ─────────────────────────────────────────────────────────────
# MVLIMS (true to Fortran: upper = max(upper, lower))
# ─────────────────────────────────────────────────────────────
@inline function _mvlins(a::T, b::T, infin::Int) where {T<:Real}
    lo = zero(T); up = one(T)
    if infin ≥ 0
        if infin != 0; lo = _cdf(a); end
        if infin != 1; up = _cdf(b); end
    end
    up = max(up, lo)
    return lo, up
end

# ─────────────────────────────────────────────────────────────
# Richtmyer roots (prime roots)
# returns [sqrt(p1), …, sqrt(p_{n-1})] typed as T
# ─────────────────────────────────────────────────────────────
@inline function richtmyer_roots(::Type{T}, n::Int) where {T<:Real}
    ub = max(n-1, Int(floor(5 * n * log(n + 1) / 4)))   # NOTE: 5*n*log(...)
    ps = Primes.primes(1, ub)
    return sqrt.(T.(float.(ps[1:n-1])))
end

# ─────────────────────────────────────────────────────────────
# O(N) Swaps in Packaging (MVSWAP)
# ─────────────────────────────────────────────────────────────
@inline function _mvswap!(p::Int, q::Int,
                          A::AbstractVector, B::AbstractVector,
                          DL::AbstractVector, INFI::AbstractVector{Int},
                          COV::AbstractVector, n::Int)
    p==q && return
    A[p],A[q] = A[q],A[p];  B[p],B[q] = B[q],B[p]
    DL[p],DL[q] = DL[q],DL[p];  INFI[p],INFI[q] = INFI[q],INFI[p]
    @inbounds begin
        COV[_lin(p,p)], COV[_lin(q,q)] = COV[_lin(q,q)], COV[_lin(p,p)]
        for j in 1:p-1
            COV[_lin(p,j)], COV[_lin(q,j)] = COV[_lin(q,j)], COV[_lin(p,j)]
        end
        for i in p+1:q-1
            COV[_lin(i,p)], COV[_lin(q,i)] = COV[_lin(q,i)], COV[_lin(i,p)]
        end
        for i in q+1:n
            COV[_lin(i,p)], COV[_lin(i,q)] = COV[_lin(i,q)], COV[_lin(i,p)]
        end
    end
end

# ─────────────────────────────────────────────────────────────
# MVSORT Core (SoA + Vector Packing)
# ─────────────────────────────────────────────────────────────
function mvsort!(A::AbstractVector{T}, B::AbstractVector{T},
                 DL::AbstractVector{T}, INFI::AbstractVector{Int},
                 COV::AbstractVector{T}, n::Int;
                 pivot::Bool=true, eps::Real=1e-10) where {T<:Real}

    @assert length(A)==n==length(B)==length(DL)==length(INFI)
    @assert length(COV)==n*(n+1)÷2

    ND = n
    @inbounds for i in 1:n
        if INFI[i] < 0; ND -= 1; end
    end

    if ND > 0
        for i in n:-1:(ND+1)
            if INFI[i] ≥ 0
                for j in 1:i-1
                    if INFI[j] < 0
                        _mvswap!(j, i, A, B, DL, INFI, COV, n)
                        break
                    end
                end
            end
        end
    end

    Y = zeros(T, n)
    INFORM = 0

    if ND > 0
        for i in 1:ND
            DEMIN = one(T); JMIN = i
            CVDIAG = zero(T); AMIN=zero(T); BMIN=zero(T)
            EPSI = T(eps)*T(i)
            JL = pivot ? ND : i

            @inbounds for j in i:JL
                cjj = COV[_lin(j,j)]
                if cjj > EPSI
                    SUMSQ = sqrt(cjj)
                    SUM = DL[j]
                    for k in 1:i-1
                        SUM += COV[_lin(j,k)]*Y[k]
                    end
                    AJ = (A[j]-SUM)/SUMSQ
                    BJ = (B[j]-SUM)/SUMSQ
                    D,E = _mvlins(AJ,BJ,INFI[j])
                    mass = E-D
                    if DEMIN ≥ mass
                        DEMIN = mass; JMIN = j
                        AMIN = AJ; BMIN = BJ; CVDIAG = SUMSQ
                    end
                end
            end

            if JMIN > i
                _mvswap!(i, JMIN, A, B, DL, INFI, COV, n)
            end

            if COV[_lin(i,i)] < -EPSI
                INFORM = 3
            end

            COV[_lin(i,i)] = CVDIAG

            if CVDIAG > 0
                @inbounds for ℓ in i+1:ND
                    COV[_lin(ℓ,i)] /= CVDIAG
                    for j in i+1:ℓ
                        COV[_lin(ℓ,j)] -= COV[_lin(ℓ,i)]*COV[_lin(j,i)]
                    end
                end
                if DEMIN > EPSI
                    yv = zero(T)
                    if INFI[i] != 0; yv += _pdf(AMIN); end
                    if INFI[i] != 1; yv -= _pdf(BMIN); end
                    Y[i] = yv/DEMIN
                else
                    Y[i] = (INFI[i]==0 ? BMIN : INFI[i]==1 ? AMIN : (AMIN+BMIN)/2)
                end
                @inbounds for j in 1:i
                    COV[_lin(i,j)] /= CVDIAG
                end
                A[i]  /= CVDIAG
                B[i]  /= CVDIAG
                DL[i] /= CVDIAG
            else
                @inbounds for ℓ in i+1:ND
                    COV[_lin(ℓ,i)] = zero(T)
                end
                @inbounds for j in i-1:-1:1
                    cij = COV[_lin(i,j)]
                    if abs(cij) > EPSI
                        A[i]  /= cij; B[i]  /= cij; DL[i] /= cij
                        if cij < 0
                            A[i],B[i] = B[i],A[i]
                            if INFI[i] != 2; INFI[i] = 1 - INFI[i]; end
                        end
                        for k in 1:j
                            COV[_lin(i,k)] /= cij
                        end
                        @inbounds for ℓ in j+1:i-1
                            if COV[_lin(ℓ,j)] > zero(T)
                                for k in i-1:-1:ℓ
                                    _mvswap!(k, k+1, A, B, DL, INFI, COV, n)
                                end
                                break
                            end
                        end
                        break
                    else
                        COV[_lin(i,j)] = zero(T)
                    end
                end
                Y[i] = zero(T)
            end
        end
    end
    return (ND, Y, INFORM)
end

# ─────────────────────────────────────────────────────────────
# Preparation from Σ (covariance or correlation)
# ─────────────────────────────────────────────────────────────
function mvprep(Σ::AbstractMatrix{T},
                a::AbstractVector{T}, b::AbstractVector{T};
                δ::AbstractVector{T}=zeros(T, size(Σ,1)),
                assume_correlation::Bool=false,
                pivot::Bool=true, eps::Real=1e-10) where {T<:Real}

    n = size(Σ,1)
    @assert size(Σ)==(n,n)
    @assert length(a)==n==length(b)==length(δ)

    # Normalize to correlation and scale limits
    if assume_correlation
        R  = Σ
        A  = copy(a);  B = copy(b);  DL = copy(δ)
    else
        σ = sqrt.(LinearAlgebra.diag(Σ))
        R = Σ ./ (σ*σ')
        A = a ./ σ; B = b ./ σ; DL = δ ./ σ
    end

    # Packed with diag=1
    COV = Vector{T}(undef, n*(n+1)÷2)
    @inbounds begin
        for i in 1:n
            for j in 1:i-1
                COV[_lin(i,j)] = R[i,j]
            end
            COV[_lin(i,i)] = one(T)
        end
    end

    # INFI from (A,B)
    INFI = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        if  isinf(A[i]) &&  isinf(B[i])
            INFI[i] = -1
        elseif isinf(A[i]) && !isinf(B[i])
            INFI[i] = 0
        elseif !isinf(A[i]) &&  isinf(B[i])
            INFI[i] = 1
        else
            INFI[i] = 2
        end
    end

    ND, Y, INFORM = mvsort!(A, B, DL, INFI, COV, n; pivot=pivot, eps=eps)
    PM = PackedMatrix{T}(COV, n)  # wrapper over the COV vector

    return (nd=ND, A=A, B=B, DL=DL, INFI=INFI, COV=PM, Y=Y, inform=INFORM)
end

# Map [0,1) → (nextfloat(0), prevfloat(1)) to avoid edges
@inline function _open01(u::T) where {T<:Real}
    lo = nextfloat(zero(T))
    hi = prevfloat(one(T))
    return lo + (hi - lo)*u
end

# Normal safe quantile at (0,1), avoids ±Inf for p=0 or p=1
@inline function _safe_qf_std(p::T) where {T<:Real}
    lo = nextfloat(zero(T))
    hi = prevfloat(one(T))
    return Distributions.quantile(Distributions.Normal(zero(T), one(T)),
                                  clamp(p, lo, hi))
end

function _qmc_integrate(nd::Int, f_builder::Function; ν::Real, maxpts::Int,
                        abseps::Real, releps::Real, T::Type{<:Real},
                        rng=Random.default_rng(), antithetic::Bool=true)

    dim = nd + (ν > 0 ? 1 : 0)
    if dim == 0
        return (value=one(T), error=zero(T), inform=0)
    end

    # Richtmyer: we ask for dim+1 and use the first dim
    α = richtmyer_roots(T, dim + 1)

    # fixed shifts (12) — good compromise for error estimation
    n_shifts = 12
    per_iter_calls = antithetic ? 2 : 1
    m = max(1, maxpts ÷ (n_shifts * per_iter_calls))

    vals  = Vector{T}(undef, n_shifts)
    shift = Vector{T}(undef, dim)
    u     = Vector{T}(undef, dim)
    v     = Vector{T}(undef, dim)

    # Construct the integrand with preassigned y (closure without internal assignments)
    ybuf = Vector{T}(undef, nd)
    f = f_builder(ybuf)  #f::Function that accepts a vector u::Vector{T}

    @inbounds for s in 1:n_shifts
        for k in 1:dim
            shift[k] = rand(rng, T)
        end

        acc = zero(T)
        for r in 1:m
            # u = frac(r*α + shift) mapped to (0,1)
            @inbounds for k in 1:dim
                tmp = r*α[k] + shift[k]
                tmp -= floor(tmp)      # frac ∈ [0,1)
                u[k] = _open01(T(tmp)) # (0,1)
            end
            acc += f(u)

            if antithetic
                # Antithetical in the first nd coordinates (those passing through Φ⁻¹)
                @inbounds for k in 1:nd
                    v[k] = _open01(one(T) - u[k])
                end
                if dim > nd
                    v[dim] = u[dim]    # not antithetical in the χ² coord
                end
                acc += f(v)
            end
        end

        vals[s] = acc / T(m * per_iter_calls)
    end

    μ = mean(vals)
    σ = (n_shifts > 1) ? std(vals) : zero(T)
    tcrit = (n_shifts > 1) ? T(quantile(TDist(n_shifts-1), 0.995)) : zero(T)
    err = tcrit * σ / sqrt(T(n_shifts))

    tol = max(T(abseps), abs(T(releps)*μ))
    inform = (err <= tol) ? 0 : 1
    return (value=μ, error=err, inform=inform)
end

@inline function _integrand!(y::AbstractVector{T}, u::AbstractVector{T}, nd::Int,
                             A::AbstractVector{T}, B::AbstractVector{T}, DL::AbstractVector{T},
                             INFI::AbstractVector{Int}, COV::PackedMatrix{T}, ν::Real) where {T<:Real}

    R = (ν > 0) ? _scale_t(ν, u[end]) : one(T)
    val = one(T)

    @inbounds for k in 1:nd
        s = DL[k]
        @inbounds for j in 1:k-1
            s += _get_cov(COV, k, j) * y[j]
        end

        ak = R*A[k] - s
        bk = R*B[k] - s
        D,E = _mvlins(ak, bk, INFI[k])
        mass = E - D
        if mass <= 0
            return zero(T)
        end
        val *= mass

        p = D + u[k]*mass
        y[k] = _safe_qf_std(p)
    end

    return val
end




# ─────────────────────────────────────────────────────────────
# API principal
# ─────────────────────────────────────────────────────────────
"""
    mvtcdf(Σ, a, b; ν=0, δ=zeros, maxpts=1000n, abseps=1e-6, releps=1e-6,
           assume_correlation=false, pivot=true, rng=Random.default_rng())

Probabilidad MVN/MVT (no central) en `[a,b]` usando
MVSORT + QMC (Richtmyer + shifts). 

Retorna `(value, error, inform)`:
- `inform = 0` éxito; `1` tol. no alcanzada con `maxpts`;
  `2` dimensión inválida; `3` matriz no-PSD.
"""
# ─────────────────────────────────────────────────────────────
# API principal
# ─────────────────────────────────────────────────────────────
"""
    mvtcdf(Σ, a, b; ν=0, δ=zeros, maxpts=1000n, abseps=1e-6, releps=1e-6,
           assume_correlation=false, pivot=true, antithetic=false,
           rng=Random.default_rng())

Probabilidad MVN/MVT (no central) en `[a,b]` usando
MVSORT + QMC (Richtmyer + shifts).

Retorna `(value, error, inform)`:
- `inform = 0` éxito; `1` tol. no alcanzada con `maxpts`;
  `2` dimensión inválida; `3` matriz no-PSD.
"""
function mvtcdf(Σ::AbstractMatrix{T},
                a::AbstractVector{T},
                b::AbstractVector{T};
                ν::Real=0,
                δ::AbstractVector{T}=zeros(T, size(Σ,1)),
                maxpts::Int=1000*length(a),
                abseps::Real=1e-6,
                releps::Real=1e-6,
                assume_correlation::Bool=false,
                pivot::Bool=true,
                antithetic::Bool=false,
                rng=Random.default_rng()) where {T<:Real}

    n = size(Σ,1)
    if n < 1 || n > 1000
        return (value=T(0), error=one(T), inform=2)
    end
    @assert size(Σ) == (n,n)
    @assert length(a) == n == length(b) == length(δ)

    res = mvprep(Σ, a, b; δ=δ, assume_correlation=assume_correlation,
                 pivot=pivot, eps=1e-10)
    if res.inform == 3
        return (value=T(0), error=one(T), inform=3)
    end
    if res.nd == 0
        return (value=one(T), error=zero(T), inform=0)
    end
    if res.nd == 1 && (ν < 1 || iszero(res.DL[1]))
        v = one(T)
        if res.INFI[1] != 1
            v = (ν < 1) ? _cdf(res.B[1] - res.DL[1]) :
                          T(cdf(TDist(ν), res.B[1] - res.DL[1]))
        end
        if res.INFI[1] != 0
            v -= (ν < 1) ? _cdf(res.A[1] - res.DL[1]) :
                           T(cdf(TDist(ν), res.A[1] - res.DL[1]))
        end
        v = max(zero(T), v)
        return (value=v, error=T(2e-16), inform=0)
    end

    nd  = res.nd
    A   = res.A;  B = res.B;  DL = res.DL
    INFI= res.INFI; COV = res.COV

    f_builder = (ybuf)->(u)->_integrand!(ybuf, u, nd, A, B, DL, INFI, COV, ν)

    return _qmc_integrate(nd, f_builder; ν=ν, maxpts=maxpts,
                          abseps=abseps, releps=releps, T=T,
                          rng=rng, antithetic=antithetic)
end