# ─────────────────────────────────────────────────────────────
# Multivariate rectangular probabilities: MVN / MVT
#
# Pure Julia implementation of randomized Genz-style rank-1 lattice
# quasi-Monte Carlo scheme. The core is specialized for the two relevant
# cases:
#   - ν <= 0 : multivariate Gaussian
#   - ν > 0  : multivariate Student t
#
# Main numerical changes relative to the previous implementation:
#   1. MVN integration dimension is nd - 1, not nd.
#   2. MVT integration dimension is nd, not nd + 1.
#   3. No closure-based integrand dispatch inside the QMC loop.
#   4. QMC generation and integrand evaluation are fused in one pass.
#   5. MVN and MVT have separate hot paths.
#   6. Independent MVN rectangles are evaluated exactly.
#   7. MVN uses a cached CBC rank-1 lattice with tent transformation by default.
#   8. MVT uses a batched reduced-dimension Genz transform with folded
#      conditional Normal coordinates.
# ─────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────
# Standard Normal helpers
# ─────────────────────────────────────────────────────────────
@inline _Z(T::Type{<:Real}) = Distributions.Normal(zero(T), one(T))
@inline _cdf(x::T) where {T<:Real} = Distributions.cdf(_Z(T), x)
@inline _pdf(x::T) where {T<:Real} = Distributions.pdf(_Z(T), x)
@inline _qf_std(p::T) where {T<:Real} = Distributions.quantile(_Z(T), p)

# χ² → scale for t: R = sqrt(χ²_ν / ν).
# The hot MVT path passes the Chisq distribution and invsqrtν explicitly to
# avoid constructing Chisq(ν) at every QMC point.
# Exact radial transform for ν = 4.
#
# If X ~ χ²₄, then Y = X/2 ~ Gamma(2, 1), with
#
#     F_Y(y) = 1 - exp(-y)(1+y).
#
# Solving F_Y(y) = u gives
#
#     y = -1 - W₋₁(-(1-u)/e),
#
# and therefore
#
#     sqrt(X/4) = sqrt(y/2).
#
# Near u = 0, evaluate 1 + W₋₁ close to its branch point
# directly with lambertwbp to avoid catastrophic cancellation.
@inline function _scale_t_nu4(χ::Distributions.Chisq, invsqrtν::T, u::T,) where {T<:AbstractFloat}

    # Extremely small probabilities are essentially impossible in the
    # randomized QMC loop, but using the generic χ² quantile here preserves
    # full robustness when u approaches the smallest representable number.
    if u < sqrt(floatmin(T))
        return sqrt(T(Distributions.quantile(χ, u))) * invsqrtν
    end

    inv_e = exp(-one(T))

    if u <= T(0.5)
        # lambertwbp(z, -1) = 1 + W₋₁(-1/e + z),
        # and -(1-u)/e = -1/e + u/e.
        v = LambertW.lambertwbp(u * inv_e, -1)
        return sqrt(-T(v) / T(2))
    else
        w = LambertW.lambertw(
            -(one(T) - u) * inv_e,
            -1,
        )
        return sqrt((-one(T) - T(w)) / T(2))
    end
end

@inline function _scale_t(ν::Real, u::T) where {T<:AbstractFloat}
    @assert ν > 0

    if ν == 1
        # χ²₁ quantile:
        # sqrt(Q(u)) = sqrt(2) * erfinv(u).
        return sqrt(T(2)) * T(SpecialFunctions.erfinv(u))
    elseif ν == 2
        # χ²₂ ~ Exponential(scale=2):
        # sqrt(Q(u)/2) = sqrt(-log(1-u)).
        return sqrt(-log1p(-u))
    elseif ν == 4
        χ = Distributions.Chisq(ν)
        return _scale_t_nu4(χ, inv(sqrt(T(ν))), u,)
    end

    return sqrt(T(Distributions.quantile(Distributions.Chisq(ν), u)) / T(ν))
end

@inline function _scale_t(χ::Distributions.Chisq, invsqrtν::T, u::T,) where {T<:AbstractFloat}

    ν = Distributions.dof(χ)

    if ν == 1
        return sqrt(T(2)) * T(SpecialFunctions.erfinv(u))
    elseif ν == 2
        return sqrt(-log1p(-u))
    elseif ν == 4
        return _scale_t_nu4(χ, invsqrtν, u)
    end

    return sqrt(T(Distributions.quantile(χ, u))) * invsqrtν
end

# ─────────────────────────────────────────────────────────────
# Lower triangular packed storage
# ─────────────────────────────────────────────────────────────
struct PackedMatrix{T<:Real}
    data::Vector{T}
    n::Int
end

# 1-based index for the symmetric packed lower triangle.
@inline _lin(i::Int, j::Int) = i ≥ j ? (i * (i - 1)) ÷ 2 + j : (j * (j - 1)) ÷ 2 + i

@propagate_inbounds @inline function _get_cov(C::PackedMatrix, i::Int, j::Int)
    return C.data[_lin(i, j)]
end

# ─────────────────────────────────────────────────────────────
# MVLIMS: conditional Normal probability limits
#
# INFI convention:
#   -1 : (-∞, +∞)
#    0 : (-∞, b]
#    1 : [a, +∞)
#    2 : [a, b]
# ─────────────────────────────────────────────────────────────
@inline function _mvlins(a::T, b::T, infin::Int) where {T<:Real}
    lo = zero(T)
    up = one(T)

    if infin ≥ 0
        if infin != 0
            lo = _cdf(a)
        end
        if infin != 1
            up = _cdf(b)
        end
    end

    # This follows the classical Genz implementation: negative masses are
    # truncated to zero by forcing upper >= lower.
    up = max(up, lo)
    return lo, up
end

# ─────────────────────────────────────────────────────────────
# Richtmyer roots
#
# For q integration coordinates, call richtmyer_roots(T, q + 1).
# This returns q values sqrt(p_1), …, sqrt(p_q).
# ─────────────────────────────────────────────────────────────
@inline function richtmyer_roots(::Type{T}, n::Int) where {T<:Real}
    n <= 1 && return T[]
    ub = max(n - 1, Int(floor(5 * n * log(n + 1) / 4)))
    ps = Primes.primes(1, ub)

    # Very defensive fallback for small/edge cases.
    while length(ps) < n - 1
        ub *= 2
        ps = Primes.primes(1, ub)
    end

    roots = Vector{T}(undef, n - 1)
    @inbounds for i in eachindex(roots)
        roots[i] = sqrt(T(ps[i]))
    end
    return roots
end

# ─────────────────────────────────────────────────────────────
# Swap variables in packed representation
# ─────────────────────────────────────────────────────────────
@inline function _mvswap!(p::Int, q::Int, A::AbstractVector, B::AbstractVector, DL::AbstractVector,
                          INFI::AbstractVector{Int}, COV::AbstractVector, n::Int)
    p == q && return

    A[p], A[q] = A[q], A[p]
    B[p], B[q] = B[q], B[p]
    DL[p], DL[q] = DL[q], DL[p]
    INFI[p], INFI[q] = INFI[q], INFI[p]

    @inbounds begin
        COV[_lin(p, p)], COV[_lin(q, q)] = COV[_lin(q, q)], COV[_lin(p, p)]

        for j in 1:p-1
            COV[_lin(p, j)], COV[_lin(q, j)] = COV[_lin(q, j)], COV[_lin(p, j)]
        end
        for i in p+1:q-1
            COV[_lin(i, p)], COV[_lin(q, i)] = COV[_lin(q, i)], COV[_lin(i, p)]
        end
        for i in q+1:n
            COV[_lin(i, p)], COV[_lin(i, q)] = COV[_lin(i, q)], COV[_lin(i, p)]
        end
    end

    return nothing
end

# ─────────────────────────────────────────────────────────────
# MVSORT / variable reordering and Cholesky-like factorization
#
# This is the Genz variable sorting step. It moves fully infinite variables
# to the end, then greedily pivots by the smallest conditional probability.
# ─────────────────────────────────────────────────────────────
function mvsort!(A::AbstractVector{T}, B::AbstractVector{T}, DL::AbstractVector{T}, INFI::AbstractVector{Int},
                 COV::AbstractVector{T}, n::Int; pivot::Bool=true, eps::Real=1e-10) where {T<:Real}

    @assert length(A) == n == length(B) == length(DL) == length(INFI)
    @assert length(COV) == n * (n + 1) ÷ 2

    nd = n
    @inbounds for i in 1:n
        if INFI[i] < 0
            nd -= 1
        end
    end

    # Move non-active variables to the end.
    if nd > 0
        for i in n:-1:(nd + 1)
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

    y = zeros(T, n)
    inform = 0

    if nd > 0
        for i in 1:nd
            demin = one(T)
            jmin = i
            cvdiag = zero(T)
            amin = zero(T)
            bmin = zero(T)
            epsi = T(eps) * T(i)
            jl = pivot ? nd : i

            @inbounds for j in i:jl
                cjj = COV[_lin(j, j)]
                if cjj > epsi
                    sumsq = sqrt(cjj)
                    s = DL[j]
                    for k in 1:i-1
                        s += COV[_lin(j, k)] * y[k]
                    end

                    aj = (A[j] - s) / sumsq
                    bj = (B[j] - s) / sumsq
                    d, e = _mvlins(aj, bj, INFI[j])
                    mass = e - d

                    if demin ≥ mass
                        demin = mass
                        jmin = j
                        amin = aj
                        bmin = bj
                        cvdiag = sumsq
                    end
                end
            end

            if jmin > i
                _mvswap!(i, jmin, A, B, DL, INFI, COV, n)
            end

            if COV[_lin(i, i)] < -epsi
                inform = 3
            end

            COV[_lin(i, i)] = cvdiag

            if cvdiag > 0
                @inbounds for ℓ in i+1:nd
                    COV[_lin(ℓ, i)] /= cvdiag
                    for j in i+1:ℓ
                        COV[_lin(ℓ, j)] -= COV[_lin(ℓ, i)] * COV[_lin(j, i)]
                    end
                end

                if demin > epsi
                    yy = zero(T)
                    if INFI[i] != 0
                        yy += _pdf(amin)
                    end
                    if INFI[i] != 1
                        yy -= _pdf(bmin)
                    end
                    y[i] = yy / demin
                else
                    y[i] = INFI[i] == 0 ? bmin : INFI[i] == 1 ? amin : (amin + bmin) / 2
                end

                @inbounds for j in 1:i
                    COV[_lin(i, j)] /= cvdiag
                end
                A[i] /= cvdiag
                B[i] /= cvdiag
                DL[i] /= cvdiag
            else
                @inbounds for ℓ in i+1:nd
                    COV[_lin(ℓ, i)] = zero(T)
                end

                @inbounds for j in i-1:-1:1
                    cij = COV[_lin(i, j)]
                    if abs(cij) > epsi
                        A[i] /= cij
                        B[i] /= cij
                        DL[i] /= cij

                        if cij < 0
                            A[i], B[i] = B[i], A[i]
                            if INFI[i] != 2
                                INFI[i] = 1 - INFI[i]
                            end
                        end

                        for k in 1:j
                            COV[_lin(i, k)] /= cij
                        end

                        for ℓ in j+1:i-1
                            if COV[_lin(ℓ, j)] > zero(T)
                                for k in i-1:-1:ℓ
                                    _mvswap!(k, k + 1, A, B, DL, INFI, COV, n)
                                end
                                break
                            end
                        end
                        break
                    else
                        COV[_lin(i, j)] = zero(T)
                    end
                end
                y[i] = zero(T)
            end
        end
    end

    return (nd, y, inform)
end

# ─────────────────────────────────────────────────────────────
# Preparation from covariance/correlation matrix
# ─────────────────────────────────────────────────────────────
function mvprep(Σ::AbstractMatrix{T}, a::AbstractVector{T}, b::AbstractVector{T};
                δ::AbstractVector{T}=zeros(T, size(Σ, 1)), assume_correlation::Bool=false,
                pivot::Bool=true, eps::Real=1e-10) where {T<:Real}

    n = size(Σ, 1)
    @assert size(Σ) == (n, n)
    @assert length(a) == n == length(b) == length(δ)

    A = Vector{T}(undef, n)
    B = Vector{T}(undef, n)
    DL = Vector{T}(undef, n)
    COV = Vector{T}(undef, n * (n + 1) ÷ 2)

    if assume_correlation
        @inbounds begin
            for i in 1:n
                A[i] = a[i]
                B[i] = b[i]
                DL[i] = δ[i]
                for j in 1:i-1
                    COV[_lin(i, j)] = Σ[i, j]
                end
                COV[_lin(i, i)] = one(T)
            end
        end
    else
        σ = Vector{T}(undef, n)
        @inbounds for i in 1:n
            σ[i] = sqrt(Σ[i, i])
        end

        @inbounds begin
            for i in 1:n
                A[i] = a[i] / σ[i]
                B[i] = b[i] / σ[i]
                DL[i] = δ[i] / σ[i]

                for j in 1:i-1
                    COV[_lin(i, j)] = Σ[i, j] / (σ[i] * σ[j])
                end
                COV[_lin(i, i)] = one(T)
            end
        end
    end

    INFI = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        if isinf(A[i]) && isinf(B[i])
            INFI[i] = -1
        elseif isinf(A[i]) && !isinf(B[i])
            INFI[i] = 0
        elseif !isinf(A[i]) && isinf(B[i])
            INFI[i] = 1
        else
            INFI[i] = 2
        end
    end

    nd, y, inform = mvsort!(A, B, DL, INFI, COV, n; pivot=pivot, eps=eps)
    return (nd=nd, A=A, B=B, DL=DL, INFI=INFI, COV=PackedMatrix{T}(COV, n), Y=y, inform=inform)
end

# ─────────────────────────────────────────────────────────────
# Small utilities
# ─────────────────────────────────────────────────────────────
@inline function _open01(u::T) where {T<:Real}
    lo = nextfloat(zero(T))
    hi = prevfloat(one(T))
    return lo + (hi - lo) * u
end

@inline function _safe_qf_std(p::T) where {T<:Real}
    lo = nextfloat(zero(T))
    hi = prevfloat(one(T))
    return _qf_std(clamp(p, lo, hi))
end

@inline function _rqmc_raw(r::Int, α::T, shift::T) where {T<:Real}
    u = T(r) * α + shift
    return u - floor(u)
end

@inline _rqmc_coord(r::Int, α::T, shift::T) where {T<:Real} = _open01(_rqmc_raw(r, α, shift))


# Folded/tent-transformed rank-1 lattice coordinate used by the Genz core.
# This maps frac(rα + shift) to |2u - 1|. It gives the same support (0,1),
# but usually lowers variation of the transformed integrand for MVN rectangles.
@inline function _rqmc_folded_coord(r::Int, α::T, shift::T) where {T<:Real}
    u = _rqmc_raw(r, α, shift)
    return _open01(abs(T(2) * u - one(T)))
end

@inline _rqmc_coord_antithetic(r::Int, α::T, shift::T) where {T<:Real} = _open01(one(T) - _rqmc_raw(r, α, shift))

@inline function _qmc_error(vals::AbstractVector{T}) where {T<:Real}
    n = length(vals)
    n <= 1 && return zero(T)
    # Classical Genz-Bretz style: approximately three standard errors.
    return T(3) * Statistics.std(vals) / sqrt(T(n))
end

@inline function _qmc_inform(value::T, error::T, abseps::Real, releps::Real) where {T<:Real}
    tol = max(T(abseps), abs(T(releps) * value))
    return error <= tol ? 0 : 1
end

function _is_diagonal(Σ::AbstractMatrix{T}, n::Int) where {T<:Real}
    @inbounds for j in 1:n
        for i in j+1:n
            if !iszero(Σ[i, j]) || !iszero(Σ[j, i])
                return false
            end
        end
    end
    return true
end

function _mvn_independent_cdf(Σ::AbstractMatrix{T}, a::AbstractVector{T}, b::AbstractVector{T},
                              δ::AbstractVector{T}; assume_correlation::Bool=false) where {T<:Real}
    n = length(a)
    value = one(T)

    @inbounds for i in 1:n
        σ = assume_correlation ? one(T) : sqrt(Σ[i, i])
        lo = (a[i] - δ[i]) / σ
        up = (b[i] - δ[i]) / σ
        mass = _cdf(up) - _cdf(lo)
        mass = clamp(mass, zero(T), one(T))
        value *= mass
    end

    return value
end

# ─────────────────────────────────────────────────────────────
# Hot integrands: reduced-dimension Genz transforms
#
# MVN uses q = nd - 1 quasi-random coordinates.
# MVT uses q = nd coordinates: 1 for the radial χ² scale and nd - 1 for
# the conditional normal transform.
# ─────────────────────────────────────────────────────────────
function _mvn_integrand_rqmc!(y::AbstractVector{T}, r::Int, α::AbstractVector{T}, shift::AbstractVector{T},
                              nd::Int, A::AbstractVector{T}, B::AbstractVector{T}, DL::AbstractVector{T},
                              INFI::AbstractVector{Int}, COV::PackedMatrix{T}, antithetic::Bool) where {T<:Real}
    value = one(T)

    @inbounds begin
        # First conditional probability is exact. Its random coordinate is
        # only needed to generate y[1] for the next conditional levels.
        s = DL[1]
        d, e = _mvlins(A[1] - s, B[1] - s, INFI[1])
        mass = e - d
        mass <= 0 && return zero(T)
        value *= mass

        if nd > 1
            u = antithetic ? _rqmc_coord_antithetic(r, α[1], shift[1]) :
                             _rqmc_folded_coord(r, α[1], shift[1])
            y[1] = _safe_qf_std(d + u * mass)
        end

        for k in 2:nd
            s = DL[k]
            for j in 1:k-1
                s += _get_cov(COV, k, j) * y[j]
            end

            d, e = _mvlins(A[k] - s, B[k] - s, INFI[k])
            mass = e - d
            mass <= 0 && return zero(T)
            value *= mass

            if k < nd
                u = antithetic ? _rqmc_coord_antithetic(r, α[k], shift[k]) :
                                 _rqmc_folded_coord(r, α[k], shift[k])
                y[k] = _safe_qf_std(d + u * mass)
            end
        end
    end

    return value
end

function _mvt_integrand_rqmc!(y::AbstractVector{T}, r::Int, α::AbstractVector{T}, shift::AbstractVector{T},
                              nd::Int, A::AbstractVector{T}, B::AbstractVector{T}, DL::AbstractVector{T},
                              INFI::AbstractVector{Int}, COV::PackedMatrix{T}, χ::Distributions.Chisq,
                              invsqrtν::T, antithetic::Bool) where {T<:Real}
    value = one(T)

    @inbounds begin
        # First QMC coordinate is the common t scale. The same tent
        # periodization used by the conditional Normal coordinates is applied
        # to the radial coordinate for the CBC lattice.
        uχ = _rqmc_folded_coord(r, α[1], shift[1])
        R = _scale_t(χ, invsqrtν, uχ)

        s = DL[1]
        d, e = _mvlins(R * A[1] - s, R * B[1] - s, INFI[1])
        mass = e - d
        mass <= 0 && return zero(T)
        value *= mass

        if nd > 1
            # Conditional Normal coordinates use the same folded/tent transform
            # that stabilized the MVN core.
            u = antithetic ? _rqmc_coord_antithetic(r, α[2], shift[2]) :
                             _rqmc_folded_coord(r, α[2], shift[2])
            y[1] = _safe_qf_std(d + u * mass)
        end

        for k in 2:nd
            s = DL[k]
            for j in 1:k-1
                s += _get_cov(COV, k, j) * y[j]
            end

            d, e = _mvlins(R * A[k] - s, R * B[k] - s, INFI[k])
            mass = e - d
            mass <= 0 && return zero(T)
            value *= mass

            if k < nd
                u = antithetic ? _rqmc_coord_antithetic(r, α[k + 1], shift[k + 1]) :
                                 _rqmc_folded_coord(r, α[k + 1], shift[k + 1])
                y[k] = _safe_qf_std(d + u * mass)
            end
        end
    end

    return value
end

# ─────────────────────────────────────────────────────────────
# Randomized rank-1 lattice QMC drivers
# ─────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────
# Batched MVN randomized-lattice QMC driver
#
# The scalar Genz transform is memory-light but pays substantial overhead per
# point. This batch version evaluates many QMC points for one random shift in a
# single pass over each conditional level. It keeps the previous low-allocation
# spirit by using a bounded batch instead of storing the full nper × (nd - 1)
# matrix used by fully vectorized implementations.
# ─────────────────────────────────────────────────────────────
@inline function _default_mvn_batchsize(nper::Int, nd::Int)
    # A power-of-two batch is cache friendly. For large dimensions this keeps
    # memory moderate: 4096 × 19 Float64 ≈ 608 KiB for nd = 20.
    nd <= 3 && return min(nper, 8192)
    nd <= 10 && return min(nper, 4096)
    return min(nper, 2048)
end

function _mvn_integrand_rqmc_batch!(pv::AbstractVector{T},
                                    y::AbstractMatrix{T},
                                    work::AbstractVector{T},
                                    r0::Int,
                                    nb::Int,
                                    α::AbstractVector{T},
                                    shift::AbstractVector{T},
                                    nd::Int,
                                    A::AbstractVector{T}, B::AbstractVector{T},
                                    DL::AbstractVector{T},
                                    INFI::AbstractVector{Int},
                                    COV::PackedMatrix{T},
                                    antithetic::Bool) where {T<:Real}
    @inbounds begin
        # First conditional probability is exact. The first QMC coordinate is
        # used only to generate y₁ for the following conditional levels.
        d, e = _mvlins(A[1] - DL[1], B[1] - DL[1], INFI[1])
        mass = e - d
        mass <= 0 && return zero(T)

        @simd for i in 1:nb
            pv[i] = mass
        end

        if nd > 1
            a1 = α[1]
            sh1 = shift[1]
            @simd for i in 1:nb
                r = r0 + i - 1
                u = antithetic ? _rqmc_coord_antithetic(r, a1, sh1) :
                                 _rqmc_folded_coord(r, a1, sh1)
                y[i, 1] = _safe_qf_std(d + u * mass)
            end
        end

        for k in 2:nd
            # work[:] = DL[k] + y[:, 1:k-1] * COV[k, 1:k-1]
            @simd for i in 1:nb
                work[i] = DL[k]
            end

            row0 = _lin(k, 1)
            for j in 1:k-1
                ckj = COV.data[row0 + j - 1]
                @simd for i in 1:nb
                    work[i] += ckj * y[i, j]
                end
            end

            if k < nd
                ak = α[k]
                shk = shift[k]
                for i in 1:nb
                    if iszero(pv[i])
                        y[i, k] = zero(T)
                        continue
                    end

                    d, e = _mvlins(A[k] - work[i], B[k] - work[i], INFI[k])
                    mass = e - d
                    if mass <= 0
                        pv[i] = zero(T)
                        y[i, k] = zero(T)
                    else
                        pv[i] *= mass
                        r = r0 + i - 1
                        u = antithetic ? _rqmc_coord_antithetic(r, ak, shk) :
                                         _rqmc_folded_coord(r, ak, shk)
                        y[i, k] = _safe_qf_std(d + u * mass)
                    end
                end
            else
                for i in 1:nb
                    if iszero(pv[i])
                        continue
                    end

                    d, e = _mvlins(A[k] - work[i], B[k] - work[i], INFI[k])
                    mass = e - d
                    pv[i] = mass <= 0 ? zero(T) : pv[i] * mass
                end
            end
        end

        acc = zero(T)
        @simd for i in 1:nb
            acc += pv[i]
        end
        return acc
    end
end

function _rqmc_integrate_mvn(nd::Int, A::AbstractVector{T}, B::AbstractVector{T}, DL::AbstractVector{T},
                             INFI::AbstractVector{Int}, COV::PackedMatrix{T}; maxpts::Int, abseps::Real,
                             releps::Real, rng=Random.default_rng(), antithetic::Bool=false, nshifts::Int=10,
                             batchsize::Int=0) where {T<:Real}
    qdim = nd - 1
    qdim <= 0 && return (value=one(T), error=zero(T), inform=0)

    per_point = antithetic ? 2 : 1
    requested_nper = max(1, maxpts ÷ (nshifts * per_point),)

    α, nper = _qmc_lattice(T, qdim, requested_nper,)

    shift = Vector{T}(undef, qdim)
    vals = Vector{T}(undef, nshifts)

    bsz = batchsize > 0 ? min(batchsize, nper) : _default_mvn_batchsize(nper, nd)

    # y stores only the generated conditional Normal quantiles y₁,…,y_{nd-1}.
    # The last conditional level contributes only a probability mass.
    y = Matrix{T}(undef, bsz, qdim)
    pv = Vector{T}(undef, bsz)
    work = Vector{T}(undef, bsz)

    @inbounds for sidx in 1:nshifts
        for k in 1:qdim
            shift[k] = rand(rng, T)
        end

        acc = zero(T)
        r0 = 1
        while r0 <= nper
            nb = min(bsz, nper - r0 + 1)
            acc += _mvn_integrand_rqmc_batch!(pv, y, work, r0, nb, α, shift,
                                              nd, A, B, DL, INFI, COV, false)
            if antithetic
                acc += _mvn_integrand_rqmc_batch!(pv, y, work, r0, nb, α, shift,
                                                  nd, A, B, DL, INFI, COV, true)
            end
            r0 += nb
        end

        vals[sidx] = acc / T(nper * per_point)
    end

    value = Statistics.mean(vals)
    error = _qmc_error(vals)
    inform = _qmc_inform(value, error, abseps, releps)
    return (value=value, error=error, inform=inform)
end

@inline function _default_mvt_batchsize(nper::Int, nd::Int)
    # MVT is more expensive per point because each point needs one χ² quantile.
    # Keep batches slightly smaller than MVN while still amortizing loop overhead.
    nd <= 3 && return min(nper, 4096)
    nd <= 10 && return min(nper, 2048)
    return min(nper, 1024)
end

function _mvt_integrand_rqmc_batch!(pv::AbstractVector{T},
                                    y::AbstractMatrix{T},
                                    work::AbstractVector{T},
                                    scale::AbstractVector{T},
                                    r0::Int,
                                    nb::Int,
                                    α::AbstractVector{T},
                                    shift::AbstractVector{T},
                                    nd::Int,
                                    A::AbstractVector{T}, B::AbstractVector{T},
                                    DL::AbstractVector{T},
                                    INFI::AbstractVector{Int},
                                    COV::PackedMatrix{T},
                                    χ::Distributions.Chisq,
                                    invsqrtν::T,
                                    antithetic::Bool) where {T<:Real}
    @inbounds begin
        # First coordinate: common t scale for each QMC point.
        aχ = α[1]
        shχ = shift[1]
        for i in 1:nb
            r = r0 + i - 1
            uχ = _rqmc_folded_coord(r, aχ, shχ)
            scale[i] = _scale_t(χ, invsqrtν, uχ)
        end

        # First conditional probability.
        for i in 1:nb
            R = scale[i]
            d, e = _mvlins(R * A[1] - DL[1], R * B[1] - DL[1], INFI[1])
            mass = e - d
            if mass <= 0
                pv[i] = zero(T)
                if nd > 1
                    y[i, 1] = zero(T)
                end
            else
                pv[i] = mass
                if nd > 1
                    u = antithetic ? _rqmc_coord_antithetic(r0 + i - 1, α[2], shift[2]) :
                                     _rqmc_folded_coord(r0 + i - 1, α[2], shift[2])
                    y[i, 1] = _safe_qf_std(d + u * mass)
                end
            end
        end

        for k in 2:nd
            @simd for i in 1:nb
                work[i] = DL[k]
            end

            row0 = _lin(k, 1)
            for j in 1:k-1
                ckj = COV.data[row0 + j - 1]
                @simd for i in 1:nb
                    work[i] += ckj * y[i, j]
                end
            end

            if k < nd
                ak = α[k + 1]
                shk = shift[k + 1]
                for i in 1:nb
                    if iszero(pv[i])
                        y[i, k] = zero(T)
                        continue
                    end

                    R = scale[i]
                    d, e = _mvlins(R * A[k] - work[i], R * B[k] - work[i], INFI[k])
                    mass = e - d

                    if mass <= 0
                        pv[i] = zero(T)
                        y[i, k] = zero(T)
                    else
                        pv[i] *= mass
                        r = r0 + i - 1
                        u = antithetic ? _rqmc_coord_antithetic(r, ak, shk) :
                                         _rqmc_folded_coord(r, ak, shk)
                        y[i, k] = _safe_qf_std(d + u * mass)
                    end
                end
            else
                for i in 1:nb
                    iszero(pv[i]) && continue
                    R = scale[i]
                    d, e = _mvlins(R * A[k] - work[i], R * B[k] - work[i], INFI[k])
                    mass = e - d
                    pv[i] = mass <= 0 ? zero(T) : pv[i] * mass
                end
            end
        end

        acc = zero(T)
        @simd for i in 1:nb
            acc += pv[i]
        end
        return acc
    end
end

function _rqmc_integrate_mvt(nd::Int,
                             A::AbstractVector{T}, B::AbstractVector{T}, DL::AbstractVector{T},
                             INFI::AbstractVector{Int}, COV::PackedMatrix{T},
                             ν::Real;
                             maxpts::Int,
                             abseps::Real,
                             releps::Real,
                             rng=Random.default_rng(),
                             antithetic::Bool=false,
                             nshifts::Int=8,
                             batchsize::Int=0) where {T<:Real}
    qdim = nd
    qdim <= 0 && return (value=one(T), error=zero(T), inform=0)

    per_point = antithetic ? 2 : 1
    requested_nper = max(1, maxpts ÷ (nshifts * per_point),)

    α, nper = _qmc_lattice(T, qdim, requested_nper,)

    shift = Vector{T}(undef, qdim)
    vals = Vector{T}(undef, nshifts)

    bsz = batchsize > 0 ? min(batchsize, nper) : _default_mvt_batchsize(nper, nd)

    # Only y₁,…,y_{nd-1} are needed. The last level contributes only a mass.
    ycols = max(nd - 1, 1)
    y = Matrix{T}(undef, bsz, ycols)
    pv = Vector{T}(undef, bsz)
    work = Vector{T}(undef, bsz)
    scale = Vector{T}(undef, bsz)

    χ = Distributions.Chisq(ν)
    invsqrtν = inv(sqrt(T(ν)))

    @inbounds for sidx in 1:nshifts
        for k in 1:qdim
            shift[k] = rand(rng, T)
        end

        acc = zero(T)
        r0 = 1
        while r0 <= nper
            nb = min(bsz, nper - r0 + 1)
            acc += _mvt_integrand_rqmc_batch!(pv, y, work, scale, r0, nb, α, shift, nd,
                                              A, B, DL, INFI, COV, χ, invsqrtν, false)
            if antithetic
                acc += _mvt_integrand_rqmc_batch!(pv, y, work, scale, r0, nb, α, shift, nd,
                                                  A, B, DL, INFI, COV, χ, invsqrtν, true)
            end
            r0 += nb
        end

        vals[sidx] = acc / T(nper * per_point)
    end

    value = Statistics.mean(vals)
    error = _qmc_error(vals)
    inform = _qmc_inform(value, error, abseps, releps)
    return (value=value, error=error, inform=inform)
end

# ─────────────────────────────────────────────────────────────
# Principal internal API
# ─────────────────────────────────────────────────────────────
"""
    mvtcdf(Σ, a, b; ν=0, δ=zeros, maxpts=1000n, abseps=1e-6,
           releps=1e-6, assume_correlation=false, pivot=true,
           antithetic=false, rng=Random.default_rng(), batchsize=0,
           nshifts=nothing)

Rectangular probability for multivariate Gaussian (`ν <= 0`) and
multivariate Student t (`ν > 0`) distributions using MVSORT plus randomized
randomized rank-1 lattice quasi-Monte Carlo with cached CBC lattice construction.

When `nshifts` is not specified, the core uses 10 randomized shifts for the
Gaussian case and 8 randomized shifts for the Student t case.

Returns a named tuple `(value, error, inform)`.

`inform` codes:
- `0`: requested tolerance reached according to the internal error estimate.
- `1`: tolerance not reached with `maxpts`.
- `2`: invalid dimension.
- `3`: covariance/correlation matrix appears non positive semidefinite.
"""
function mvtcdf(Σ::AbstractMatrix{T},
                a::AbstractVector{T},
                b::AbstractVector{T};
                ν::Real=0,
                δ::AbstractVector{T}=zeros(T, size(Σ, 1)),
                maxpts::Int=1000 * length(a),
                abseps::Real=1e-6,
                releps::Real=1e-6,
                assume_correlation::Bool=false,
                pivot::Bool=true,
                antithetic::Bool=false,
                rng=Random.default_rng(),
                batchsize::Int=0,
                nshifts::Union{Nothing,Int}=nothing) where {T<:Real}

    n = size(Σ, 1)
    if n < 1 || n > 1000
        return (value=zero(T), error=one(T), inform=2)
    end

    @assert size(Σ) == (n, n)
    @assert length(a) == n == length(b) == length(δ)

    nshifts_eff = isnothing(nshifts) ? (ν > 0 ? 8 : 10) : nshifts
    if nshifts_eff < 2
        throw(ArgumentError("nshifts must be at least 2"))
    end

    # Empty rectangle. High-level wrappers already handle this, but keeping it
    # here makes the core safer when called directly.
    @inbounds for i in 1:n
        if b[i] < a[i]
            return (value=zero(T), error=zero(T), inform=0)
        end
    end

    # Exact independent Gaussian path. This matters for both accuracy and speed.
    if ν <= 0 && _is_diagonal(Σ, n)
        value = _mvn_independent_cdf(Σ, a, b, δ; assume_correlation=assume_correlation)
        return (value=value, error=zero(T), inform=0)
    end

    res = mvprep(Σ, a, b; δ=δ, assume_correlation=assume_correlation, pivot=pivot, eps=1e-10)

    if res.inform == 3
        return (value=zero(T), error=one(T), inform=3)
    end

    if res.nd == 0
        return (value=one(T), error=zero(T), inform=0)
    end

    # One-dimensional exact path.
    if res.nd == 1 && (ν <= 0 || iszero(res.DL[1]))
        value = one(T)

        if res.INFI[1] != 1
            value = ν <= 0 ? _cdf(res.B[1] - res.DL[1]) :
                             T(Distributions.cdf(Distributions.TDist(ν), res.B[1] - res.DL[1]))
        end
        if res.INFI[1] != 0
            value -= ν <= 0 ? _cdf(res.A[1] - res.DL[1]) :
                              T(Distributions.cdf(Distributions.TDist(ν), res.A[1] - res.DL[1]))
        end

        value = clamp(value, zero(T), one(T))
        return (value=value, error=T(2e-16), inform=0)
    end

    nd = res.nd
    A = res.A
    B = res.B
    DL = res.DL
    INFI = res.INFI
    COV = res.COV

    if ν > 0
        return _rqmc_integrate_mvt(nd, A, B, DL, INFI, COV, ν;
                                   maxpts=maxpts,
                                   abseps=abseps,
                                   releps=releps,
                                   rng=rng,
                                   antithetic=antithetic,
                                   nshifts=nshifts_eff,
                                   batchsize=batchsize)
    else
        return _rqmc_integrate_mvn(nd, A, B, DL, INFI, COV;
                                   maxpts=maxpts,
                                   abseps=abseps,
                                   releps=releps,
                                   rng=rng,
                                   antithetic=antithetic,
                                   nshifts=nshifts_eff,
                                   batchsize=batchsize)
    end
end
