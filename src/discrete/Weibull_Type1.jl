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
function _weibull_type1_sum(f::F, ::Type{T}, start::Integer; atol::Real = eps(T), rtol::Real = sqrt(eps(T)), maxiter::Integer = 10_000_000,) where {F,T<:AbstractFloat}

    s = zero(T)
    c = zero(T)

    atol_T = T(atol)
    rtol_T = T(rtol)

    for k in start:maxiter
        term = T(f(k))

        if term == zero(T)
            return s + c
        end

        t = s + term

        if abs(s) >= abs(term)
            c += (s - t) + term
        else
            c += (term - t) + s
        end

        s = t
        total = s + c

        if k > start && abs(term) <= max(atol_T, rtol_T * abs(total))
            return total
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
    return _weibull_type1_sum(k -> exp(logq * (T(k)^β)), T, 1; kwargs...)
end

function Statistics.var(d::Weibull_Type1; kwargs...)
    q, β = float(d.q), float(d.β)
    T = typeof(q + β)
    if β == one(T)
        return q / (one(T) - q)^2
    end
    logq = log(q)
    μ = Statistics.mean(d; kwargs...)
    m2 = _weibull_type1_sum(k -> (T(2) * T(k) - one(T)) * exp(logq * (T(k)^β)), T, 1; kwargs...)
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
    return _weibull_type1_sum(
        k -> begin
            lp = Distributions.logpdf(d, k)
            isfinite(lp) ? -exp(lp) * lp : zero(T)
        end,
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
    k = floor(x)
    logq = log(q)
    a = logq * k^β
    b = logq * ((k + 1)^β - k^β)
    # log(q^(k^β) - q^((k+1)^β))
    # = k^β log(q) + log(1 - exp(((k+1)^β - k^β)log(q))).
    return a + (b < -log(2) ? log1p(-exp(b)) : log(-expm1(b)))
end

function Distributions.quantile(d::Weibull_Type1, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    p == 0 && return 0
    p == 1 && return oftype(float(d.q), Inf)
    q, β = d.q, d.β
    return max(0, ceil(Int, (log1p(-p) / log(q))^(1 / β) - 1))
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Weibull_Type1)
    q, β = d.q, d.β
    return floor(Int, (log1p(-rand(rng)) / log(q))^(1 / β))
end