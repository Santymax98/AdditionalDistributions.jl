# ─────────────────────────────────────────────────────────────
# Fast component-by-component (CBC) rank-1 lattice construction
# ─────────────────────────────────────────────────────────────
#
# The randomized QMC integrators use a rank-1 lattice with a generating
# vector selected component-by-component. The CBC objective uses the
# Bernoulli-polynomial kernel employed by modern Genz-style implementations.
#
# Construction is accelerated with FFTs and cached by (floating type, qdim,
# prime points per random shift). The FFT work is therefore paid only on the
# first request for a given lattice.
#
# Reference:
#   Nuyens, D. and Cools, R. (2006).
#   Fast component-by-component construction, a reprise for different kernels.
#   Monte Carlo and Quasi-Monte Carlo Methods 2004, 371–385.
# ─────────────────────────────────────────────────────────────

const _CBC_LATTICE_CACHE = Dict{Tuple{DataType,Int,Int},Any}()
const _CBC_LATTICE_CACHE_LOCK = ReentrantLock()

@inline function _cbc_powermod(a::Int, e::Int, m::Int)
    result = 1
    base = mod(a, m)
    exponent = e

    while exponent > 0
        if isodd(exponent)
            result = Int(mod(Int128(result) * base, m))
        end

        exponent >>= 1

        if exponent > 0
            base = Int(mod(Int128(base) * base, m))
        end
    end

    return result
end

function _cbc_primitive_root_prime(p::Int)
    p >= 3 || throw(ArgumentError("CBC prime must be at least 3"))
    Primes.isprime(p) ||
        throw(ArgumentError("CBC modulus must be prime"))

    factors = keys(Primes.factor(p - 1))

    for g in 2:(p - 1)
        good = true

        for q in factors
            if _cbc_powermod(g, (p - 1) ÷ q, p) == 1
                good = false
                break
            end
        end

        good && return g
    end

    error("primitive root not found for prime $p")
end

function _fast_cbc_lattice_float64(qdim::Int, requested_nper::Int,)
    qdim >= 1 || throw(ArgumentError("qdim must be positive"))
    requested_nper >= 3 || throw(ArgumentError("requested_nper must be at least 3"))

    p = Int(Primes.prevprime(requested_nper))
    m = (p - 1) ÷ 2
    g = _cbc_primitive_root_prime(p)

    perm = Vector{Int}(undef, m)
    perm[1] = 1

    @inbounds for j in 2:m
        perm[j] = mod(g * perm[j - 1], p)
    end

    @inbounds for j in 1:m
        perm[j] = min(perm[j], p - perm[j])
    end

    invp = inv(Float64(p))
    kernel = Vector{Float64}(undef, m)

    @inbounds for j in 1:m
        x = perm[j] * invp
        kernel[j] = x * x - x + 1 / 6
    end

    kernel_fft = FFTW.fft(kernel)

    product_weights = ones(Float64, m)
    reordered = similar(kernel)
    z = ones(Int, qdim)

    w = 0

    for s in 2:qdim
        @inbounds for j0 in 0:(m - 1)
            reordered[j0 + 1] =
                kernel[mod(w - j0, m) + 1]
        end

        γ = s <= 3 ? 1.0 : 0.8^(s - 3)

        @inbounds @simd for j in 1:m
            product_weights[j] *= 1 + γ * reordered[j]
        end

        objective = real.(FFTW.ifft(kernel_fft .* FFTW.fft(product_weights)))

        minval = minimum(objective)
        tol = 128 * eps(Float64) * max(1.0, abs(minval))

        candidates = findall(x -> x <= minval + tol, objective)

        w = maximum(candidates) - 1
        z[s] = perm[w + 1]
    end

    α = Float64.(z) ./ p
    return α, p
end

function _cached_cbc_lattice(::Type{T}, qdim::Int, requested_nper::Int,) where {T<:AbstractFloat}
    p = Int(Primes.prevprime(requested_nper))
    key = (T, qdim, p)

    lock(_CBC_LATTICE_CACHE_LOCK)
    try
        if haskey(_CBC_LATTICE_CACHE, key)
            return _CBC_LATTICE_CACHE[key], p
        end
    finally
        unlock(_CBC_LATTICE_CACHE_LOCK)
    end

    α64, actual_p = _fast_cbc_lattice_float64(qdim, requested_nper,)

    @assert actual_p == p

    α = T.(α64)

    lock(_CBC_LATTICE_CACHE_LOCK)
    try
        cached = get!(_CBC_LATTICE_CACHE, key, α,)
        return cached, p
    finally
        unlock(_CBC_LATTICE_CACHE_LOCK)
    end
end

#     _qmc_lattice(T, qdim, requested_nper)
#
# Return `(α, nper)` for the internal randomized rank-1 lattice.
#
# For ordinary floating-point calculations and at least three requested
# points per random shift, a cached Fast-CBC lattice is used. Very small
# budgets and non-floating scalar types retain the previous Richtmyer rule
# as a compatibility fallback.
function _qmc_lattice(::Type{T}, qdim::Int, requested_nper::Int,) where {T<:AbstractFloat}
    qdim <= 0 && return T[], requested_nper

    if requested_nper < 3
        return (richtmyer_roots(T, qdim + 1), requested_nper,)
    end

    return _cached_cbc_lattice(T, qdim, requested_nper,)
end

function _qmc_lattice(::Type{T}, qdim::Int, requested_nper::Int,) where {T<:Real}
    return (richtmyer_roots(T, qdim + 1), requested_nper,)
end
