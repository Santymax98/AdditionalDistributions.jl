"""

    Weibull_Type1(q, β)

A *Type-I Discrete Weibull* distribution, introduced by Nakagawa and Osaki, is a discrete analogue of the Weibull distribution whose survival function
mimics the survival function of the continuous Weibull distribution. It has two parameters: `q`, with ``0 < q < 1``, and a shape parameter ``β > 0``.

```math
P(X = k) = q^{k^\\beta} - q^{(k+1)^\\beta}, \\qquad k = 0,1,2,\\ldots
```

The cumulative distribution function is

```math
F(k) = 1 - q^{(k+1)^\\beta}, \\qquad k = 0,1,2,\\ldots
```

When `β = 1`, `Weibull_Type1(q, β)` reduces to a geometric distribution on `\\{0,1,2,\\ldots\\}`.

```julia
Weibull_Type1()       # equivalent to Weibull_Type1(0.5, 1)
Weibull_Type1(q)      # not defined; both q and β are required
Weibull_Type1(q, β)

params(d)             # Get the parameters, i.e. (q, β)
scale(d)              # Get q
shape(d)              # Get β
```

External links

* [Discrete Weibull distribution on Wikipedia](https://en.wikipedia.org/wiki/Discrete_Weibull_distribution)
* [Vila, Nakano and Saulo (2018), Theoretical results on the discrete Weibull distribution of Nakagawa and Osaki](https://doi.org/10.1080/02331888.2018.1550645)
"""
struct Weibull_Type1{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    q::T
    β::T
    Weibull_Type1{T}(q, β) where {T<:Real} = new{T}(q, β)
end

Weibull_Type1(q::Integer, β::Integer; check_args::Bool=true) = Weibull_Type1(float(q), float(β); check_args=check_args)

function Weibull_Type1(q::Real, β::Real; check_args::Bool=true)
    @check_args Weibull_Type1 (q, zero(q) < q < one(q)) (β, β > zero(β))
    q, β = promote(q, β)
    return Weibull_Type1{typeof(q)}(q, β)
end

Weibull_Type1() = Weibull_Type1{Float64}(0.5, 1)
@distr_support Weibull_Type1 0.0 Inf

# parameters
params(d::Weibull_Type1) = (d.q, d.β)
@inline partype(d::Weibull_Type1{T}) where {T<:Real} = T
Base.eltype(::Type{Weibull_Type1{T}}) where {T} = T 

scale(d::Weibull_Type1) = d.q
shape(d::Weibull_Type1) = d.β

#statistic
@inline function _weibull_type1_logpmf(logq, β, x)
    xβ = x^β
    Δ = iszero(x) ? one(xβ) : xβ * expm1(β * log1p(inv(x)))
    a = logq * xβ
    b = logq * Δ
    return a + (b < -log(2) ? log1p(-exp(b)) : log(-expm1(b)))
end

function _weibull_type1_sum(
    f::F,
    tailf::G,
    decreasing::H,
    ::Type{T},
    start::Integer;
    atol::Real=eps(T),
    rtol::Real=sqrt(eps(T)),
    maxiter::Integer=10_000_000,
) where {F,G,H,T<:AbstractFloat}
    s = zero(T)
    c = zero(T)

    atol_T = T(atol)
    rtol_T = T(rtol)
    next_tail_check = start

    for k in start:maxiter
        term = T(f(k))

        t = s + term
        if abs(s) >= abs(term)
            c += (s - t) + term
        else
            c += (term - t) + s
        end
        s = t
        total = s + c
        tol = max(atol_T, rtol_T * abs(total))

        candidate = term == zero(T) || (k > start && abs(term) <= tol)
        if candidate && (term == zero(T) || k >= next_tail_check) && decreasing(T(k))
            # For a nonnegative decreasing tail f, the integral test gives
            # sum_{j=k+1}^∞ f(j) <= integral_k^∞ f(x) dx.
            tail, err = QuadGK.quadgk(
                tailf,
                T(k),
                T(Inf);
                atol=zero(T),
                rtol=sqrt(eps(T)),
            )
            tail_upper = abs(tail) + abs(err)
            tail_upper <= tol && return total
            next_tail_check = max(k + 1, min(maxiter, 2k))
        end

        if term == zero(T)
            throw(ArgumentError("series terms underflowed before the requested tolerance was reached"))
        end
    end

    throw(ArgumentError("series did not converge after $maxiter terms"))
end

function Statistics.mean(d::Weibull_Type1; kwargs...)
    q, β = float(d.q), float(d.β)
    T = typeof(q + β)
    if β == one(T)
        return q / (one(T) - q)
    end

    logq = log(q)
    tailf = x -> isinf(x) ? zero(T) : exp(logq * x^β)
    return _weibull_type1_sum(
        k -> tailf(T(k)),
        tailf,
        x -> true,
        T,
        1;
        kwargs...
    )
end

function Statistics.var(d::Weibull_Type1; kwargs...)
    q, β = float(d.q), float(d.β)
    T = typeof(q + β)
    if β == one(T)
        return q / (one(T) - q)^2
    end

    logq = log(q)
    c0 = -logq
    μ = Statistics.mean(d; kwargs...)

    tailf = x -> isinf(x) ? zero(T) : (T(2) * x - one(T)) * exp(logq * x^β)
    decreasing = x -> begin
        x >= one(T) || return false
        return T(2) / (T(2) * x - one(T)) <= c0 * β * x^(β - one(T))
    end

    m2 = _weibull_type1_sum(
        k -> tailf(T(k)),
        tailf,
        decreasing,
        T,
        1;
        kwargs...
    )
    σ2 = m2 - μ^2
    return σ2 < zero(T) ? zero(T) : σ2
end

Statistics.std(d::Weibull_Type1; kwargs...) = sqrt(Statistics.var(d; kwargs...))

function StatsBase.entropy(d::Weibull_Type1; kwargs...)
    q, β = float(d.q), float(d.β)
    T = typeof(q + β)
    if β == one(T)
        return -log1p(-q) - (q / (one(T) - q)) * log(q)
    end

    logq = log(q)
    c0 = -logq
    tailf = x -> begin
        isinf(x) && return zero(T)
        lp = _weibull_type1_logpmf(logq, β, x)
        return isfinite(lp) ? -exp(lp) * lp : zero(T)
    end

    decreasing = x -> begin
        lp = _weibull_type1_logpmf(logq, β, x)
        lp <= -one(T) || return false

        # For β <= 1 the continuous PMF extension is decreasing. For β > 1,
        # the condition below is sufficient for it to be decreasing from x onward.
        β <= one(T) && return true
        x > zero(T) || return false
        return x^β >= (β - one(T)) / (c0 * β)
    end

    return _weibull_type1_sum(
        k -> tailf(T(k)),
        tailf,
        decreasing,
        T,
        0;
        kwargs...
    )
end

function Distributions.cdf(d::Weibull_Type1, x::Real)
    isnan(float(x)) && return NaN
    x < 0 && return 0.0
    isinf(x) && return 1.0

    q, β = d.q, d.β
    return -expm1(log(q) * (floor(x) + 1)^β)
end

function Distributions.pdf(d::Weibull_Type1, x::Real)
    isfinite(x) || return 0.0
    x >= 0 || return 0.0
    isinteger(x) || return 0.0
    return exp(Distributions.logpdf(d, x))
end

function Distributions.logpdf(d::Weibull_Type1, x::Real)
    isfinite(x) || return -Inf
    x >= 0 || return -Inf
    isinteger(x) || return -Inf
    q, β = d.q, d.β
    return _weibull_type1_logpmf(log(q), β, floor(x))
end

function Distributions.quantile(d::Weibull_Type1, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    p == 0 && return 0
    p == 1 && return oftype(float(d.q), Inf)

    q, β = d.q, d.β
    t = (log1p(-p) / log(q))^(1 / β)
    k = max(0, floor(Int, t))

    # Floating-point rounding can place t just above or below an integer at a
    # CDF jump. Verify the candidate against the discrete quantile definition.
    while k > 0 && Distributions.cdf(d, k - 1) >= p
        k -= 1
    end
    while Distributions.cdf(d, k) < p
        k += 1
    end

    return k
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Weibull_Type1)
    q, β = d.q, d.β
    return floor(Int, (log1p(-rand(rng)) / log(q))^(1 / β))
end