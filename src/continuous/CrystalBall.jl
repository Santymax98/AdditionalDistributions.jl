"""
    CrystalBall(α, m, x̄, σ)

Crystal Ball distribution with shape parameters `α > 0`, `m > 1`, location `x̄`, and scale `σ > 0`.

The density is Gaussian to the right of `z = -α` and has a power-law left tail, where `z = (x - x̄) / σ`.

```julia
CrystalBall(α, m)        # equivalent to CrystalBall(α, m, 0, 1)
CrystalBall(α, m, x̄)     # equivalent to CrystalBall(α, m, x̄, 1)
params(d)                # returns (α, m, x̄, σ)
```
"""
struct CrystalBall{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    α::T # shape
    m::T # shape
    x̄::T # location
    σ::T # scale
    CrystalBall{T}(α, m, x̄, σ) where {T<:Real} = new{T}(α, m, x̄, σ)
end

CrystalBall(α::Integer, m::Integer, x̄::Integer, σ::Integer; check_args::Bool=true) =
    CrystalBall(float(α), float(m), float(x̄), float(σ); check_args=check_args)

function CrystalBall(α::Real, m::Real, x̄::Real, σ::Real; check_args::Bool=true)
    α, m, x̄, σ = promote(α, m, x̄, σ)
    @check_args CrystalBall (α, α > zero(α)) (m, m > one(m)) (σ, σ > zero(σ))
    return CrystalBall{typeof(α)}(α, m, x̄, σ)
end

CrystalBall(α, m) = CrystalBall(α, m, 0, 1)
CrystalBall(α, m, x̄) = CrystalBall(α, m, x̄, 1)

@distr_support CrystalBall -Inf Inf

params(d::CrystalBall) = (d.α, d.m, d.x̄, d.σ)
@inline partype(d::CrystalBall{T}) where {T<:Real} = T
Base.eltype(::Type{CrystalBall{T}}) where {T} = T

location(d::CrystalBall) = d.x̄
shape(d::CrystalBall) = (d.α, d.m)
scale(d::CrystalBall) = d.σ

A(d::CrystalBall) = (d.m / abs(d.α))^d.m * exp(-abs(d.α)^2 / 2)
B(d::CrystalBall) = d.m / abs(d.α) - abs(d.α)
C(d::CrystalBall) = d.m / (abs(d.α) * (d.m - 1)) * exp(-abs(d.α)^2 / 2)
D(d::CrystalBall) = sqrt(π / 2) * (1 + SpecialFunctions.erf(abs(d.α) / sqrt(2)))
N(d::CrystalBall) = inv(d.σ * (C(d) + D(d)))

function moments(d::CrystalBall, k::Integer)
    α, m, x̄, σ = d.α, d.m, d.x̄, d.σ
    m <= k + 1 && return NaN

    bound = x̄ - α * σ
    normal_part, _ = QuadGK.quadgk(x -> x^k * exp(-((x - x̄) / σ)^2 / 2), bound, Inf)
    power_part, _ = QuadGK.quadgk(x -> x^k * A(d) * (B(d) - (x - x̄) / σ)^(-m), -Inf, bound)

    return N(d) * (normal_part + power_part)
end

Statistics.mean(d::CrystalBall) = d.m > 2 ? moments(d, 1) : NaN
Statistics.var(d::CrystalBall) = d.m > 3 ? moments(d, 2) - moments(d, 1)^2 : NaN

function Distributions.cdf(d::CrystalBall, x::Real)
    α, m, x̄, σ = d.α, d.m, d.x̄, d.σ
    isnan(float(x)) && return NaN
    x == -Inf && return zero(partype(d))
    x == Inf && return one(partype(d))

    z = (x - x̄) / σ
    if z <= -α
        return N(d) * σ * A(d) * (B(d) - z)^(1 - m) / (m - 1)
    else
        gaussian_mass = sqrt(π / 2) * (SpecialFunctions.erf(z / sqrt(2)) + SpecialFunctions.erf(α / sqrt(2)))
        return N(d) * σ * (C(d) + gaussian_mass)
    end
end

function Distributions.pdf(d::CrystalBall, x::Real)
    α, m, x̄, σ = d.α, d.m, d.x̄, d.σ
    !isfinite(x) && return zero(partype(d))

    z = (x - x̄) / σ
    if z > -α
        return N(d) * exp(-z^2 / 2)
    else
        return N(d) * A(d) * (B(d) - z)^(-m)
    end
end

function Distributions.logpdf(d::CrystalBall, x::Real)
    α, m, x̄, σ = d.α, d.m, d.x̄, d.σ
    !isfinite(x) && return -Inf

    z = (x - x̄) / σ
    if z > -α
        return log(N(d)) - z^2 / 2
    else
        return log(N(d)) + log(A(d)) - m * log(B(d) - z)
    end
end

function Distributions.quantile(d::CrystalBall, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    p == 0 && return -Inf
    p == 1 && return Inf

    α, m, x̄, σ = d.α, d.m, d.x̄, d.σ
    tail_mass = C(d) / (C(d) + D(d))

    if p <= tail_mass
        z = B(d) - ((p * (m - 1)) / (N(d) * σ * A(d)))^(inv(1 - m))
    else
        erf_arg = (p / (N(d) * σ) - C(d)) / sqrt(π / 2) - SpecialFunctions.erf(α / sqrt(2))
        erf_arg = clamp(erf_arg, -one(erf_arg), one(erf_arg))
        z = sqrt(2) * SpecialFunctions.erfinv(erf_arg)
    end

    return x̄ + σ * z
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::CrystalBall)
    return Distributions.quantile(d, rand(rng))
end
