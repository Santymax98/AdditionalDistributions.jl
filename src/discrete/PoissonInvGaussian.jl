"""
    PoissonInvGaussian(m, ϕ)

Distribución Poisson-Inversa Gaussiana con media `m > 0` y parámetro
`ϕ > 0`.

```math
P(X=x) = \\frac{1}{x!}\\sqrt{\\frac{2m\\phi}{\\pi}}\\,e^{\\phi}
\\left(\\sqrt{\\frac{m\\phi}{2(1+\\phi/(2m))}}\\right)^{x-1/2}
K_{x-1/2}\\!\\left(\\sqrt{2m\\phi\\left(1+\\frac{\\phi}{2m}\\right)}\\right),
\\quad x=0,1,2,\\dots
```

```julia
PoissonInvGaussian()        # equivalent to PoissonInvGaussian(1.0, 1.0)
PoissonInvGaussian(m)       # equivalent to PoissonInvGaussian(m, 1.0)
PoissonInvGaussian(m, ϕ)

params(d)                   # Get the parameters, i.e. (m, ϕ)
```

Referencia: Shaban, S. A. (1981). Computation of the Poisson-inverse
Gaussian distribution. Comm. Stat. - Theory and Methods, 10(14), 1389-1399.
"""
struct PoissonInvGaussian{T<:Real} <: DiscreteUnivariateDistribution
    m::T
    ϕ::T
    PoissonInvGaussian{T}(m, ϕ) where {T<:Real} = new{T}(m, ϕ)
end

PoissonInvGaussian(m::Integer, ϕ::Integer; check_args::Bool=true) = PoissonInvGaussian(float(m), float(ϕ); check_args=check_args)

function PoissonInvGaussian(m::Real, ϕ::Real; check_args::Bool=true)
    @check_args PoissonInvGaussian (m, m > zero(m)) (ϕ, ϕ > zero(ϕ))
    m, ϕ = promote(m, ϕ)
    return PoissonInvGaussian{typeof(m)}(m, ϕ)
end

PoissonInvGaussian() = PoissonInvGaussian{Float64}(1.0, 1.0)
PoissonInvGaussian(m::Real) = PoissonInvGaussian(m, 1.0)

@distr_support PoissonInvGaussian 0 Inf

params(d::PoissonInvGaussian) = (d.m, d.ϕ)
@inline partype(d::PoissonInvGaussian{T}) where {T<:Real} = T
Base.eltype(::Type{PoissonInvGaussian{T}}) where {T} = T

Statistics.mean(d::PoissonInvGaussian) = d.m

function Statistics.var(d::PoissonInvGaussian)
    m, ϕ = params(d)
    return m + m^2 / ϕ
end

Statistics.std(d::PoissonInvGaussian) = sqrt(var(d))
Statistics.median(d::PoissonInvGaussian) = quantile(d, 0.5)

function StatsBase.skewness(d::PoissonInvGaussian)
    m, ϕ = params(d)
    κ2 = m + m^2 / ϕ
    κ3 = m + 3m^2 / ϕ + 3m^3 / ϕ^2
    return κ3 / κ2^(3 / 2)
end

function StatsBase.kurtosis(d::PoissonInvGaussian)
    m, ϕ = params(d)
    κ2 = m + m^2 / ϕ
    κ4 = m + 7m^2 / ϕ + 18m^3 / ϕ^2 + 15m^4 / ϕ^3
    return κ4 / κ2^2
end

function Distributions.cdf(d::PoissonInvGaussian, x::Real)
    isnan(x) && return NaN
    x < 0 && return 0.0
    isinf(x) && return 1.0

    n = floor(Int, x)
    m, ϕ = params(d)

    a = inv(2 * (1 + ϕ / (2m)))
    b = m * ϕ * a
    s0 = sqrt(1 + 2m / ϕ)

    p0 = exp(-2m / (1 + s0))
    n == 0 && return p0

    pkm1 = p0
    pk = p0 * sqrt(b)
    acc = p0 + pk

    for k in 1:(n - 1)
        pnext = (b / (k * (k + 1))) * pkm1 + ((2k - 1) / (k + 1)) * a * pk
        acc += pnext
        pkm1, pk = pk, pnext
    end

    return min(acc, 1.0)
end

function Distributions.pdf(d::PoissonInvGaussian, x::Real)
    insupport(d, x) || return 0.0

    n = round(Int, x)
    m, ϕ = params(d)

    a = inv(2 * (1 + ϕ / (2m)))
    b = m * ϕ * a
    s = sqrt(1 + 2m / ϕ)

    p0 = exp(-2m / (1 + s))
    n == 0 && return p0

    pkm1 = p0
    pk = p0 * sqrt(b)
    n == 1 && return pk

    for k in 1:(n - 1)
        pnext = (b / (k * (k + 1))) * pkm1 + ((2k - 1) / (k + 1)) * a * pk
        pkm1, pk = pk, pnext
    end

    return pk
end

function Distributions.logpdf(d::PoissonInvGaussian, x::Real)
    insupport(d, x) || return -Inf
    p = pdf(d, x)
    return p > 0 ? log(p) : -Inf
end

function Distributions.quantile(d::PoissonInvGaussian, p::Real)
    0 <= p <= 1 || throw(ArgumentError("p debe estar en [0, 1]"))
    p == 0 && return 0
    p == 1 && return Inf

    m, ϕ = params(d)

    a = inv(2 * (1 + ϕ / (2m)))
    b = m * ϕ * a
    s0 = sqrt(1 + 2m / ϕ)

    p0 = exp(-2m / (1 + s0))
    p <= p0 && return 0

    pkm1 = p0
    pk = p0 * sqrt(b)
    acc = p0
    k = 1

    while true
        acc += pk
        acc >= p && return k

        pnext = (b / (k * (k + 1))) * pkm1 + ((2k - 1) / (k + 1)) * a * pk
        pkm1, pk = pk, pnext
        k += 1

        k > 10^7 && throw(ErrorException("el cuantil no convergió"))
    end
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::PoissonInvGaussian)
    m, ϕ = d.m, d.ϕ
    u = rand(rng, Distributions.InverseGaussian(1, ϕ))
    return rand(rng, Distributions.Poisson(m * u))
end