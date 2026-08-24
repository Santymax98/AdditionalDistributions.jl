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

@inline function _pig_logparams(d::PoissonInvGaussian)
    m, ϕ = params(d)

    logm = log(m)
    logϕ = log(ϕ)
    log2 = log(one(m) + one(m))

    # a = m / (2m + ϕ), b = m^2 ϕ / (2m + ϕ)
    logden = LogExpFunctions.logaddexp(log2 + logm, logϕ)
    loga = logm - logden
    logb = 2logm + logϕ - logden

    s = sqrt(one(m) + 2m / ϕ)
    logp0 = -2m / (one(m) + s)

    return loga, logb, logp0
end

@inline function _pig_logstep(logpkm1, logpk, loga, logb, k::Int)
    T = typeof(logpk)
    logk = log(T(k))
    logkp1 = log(T(k + 1))

    t1 = logb - logk - logkp1 + logpkm1
    t2 = log(T(2k - 1)) - logkp1 + loga + logpk

    return LogExpFunctions.logaddexp(t1, t2)
end

function _pig_logpdf(d::PoissonInvGaussian, n::Int)
    loga, logb, logp0 = _pig_logparams(d)

    n == 0 && return logp0

    logpkm1 = logp0
    logpk = logp0 + logb / 2
    n == 1 && return logpk

    for k in 1:(n - 1)
        logpnext = _pig_logstep(logpkm1, logpk, loga, logb, k)
        logpkm1, logpk = logpk, logpnext
    end

    return logpk
end

function Distributions.cdf(d::PoissonInvGaussian, x::Real)
    isnan(x) && return NaN
    x < 0 && return 0.0
    isinf(x) && return 1.0

    n = floor(Int, x)
    loga, logb, logp0 = _pig_logparams(d)

    n == 0 && return exp(logp0)

    logpkm1 = logp0
    logpk = logp0 + logb / 2
    logacc = LogExpFunctions.logaddexp(logp0, logpk)

    for k in 1:(n - 1)
        logpnext = _pig_logstep(logpkm1, logpk, loga, logb, k)
        logacc = LogExpFunctions.logaddexp(logacc, logpnext)
        logpkm1, logpk = logpk, logpnext
    end

    return min(exp(logacc), one(logacc))
end

function Distributions.pdf(d::PoissonInvGaussian, x::Real)
    insupport(d, x) || return 0.0

    n = round(Int, x)
    return exp(_pig_logpdf(d, n))
end

function Distributions.logpdf(d::PoissonInvGaussian, x::Real)
    insupport(d, x) || return -Inf

    n = round(Int, x)
    return _pig_logpdf(d, n)
end

function Distributions.quantile(d::PoissonInvGaussian, p::Real)
    0 <= p <= 1 || throw(ArgumentError("p debe estar en [0, 1]"))
    p == 0 && return 0
    p == 1 && return Inf

    logtarget = log(p)
    loga, logb, logp0 = _pig_logparams(d)

    logtarget <= logp0 && return 0

    logpkm1 = logp0
    logpk = logp0 + logb / 2
    logacc = logp0
    k = 1

    while true
        logacc = LogExpFunctions.logaddexp(logacc, logpk)
        logacc >= logtarget && return k

        logpnext = _pig_logstep(logpkm1, logpk, loga, logb, k)
        logpkm1, logpk = logpk, logpnext
        k += 1

        k > 10^7 && throw(ErrorException("el cuantil no convergió"))
    end
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::PoissonInvGaussian)
    m, ϕ = d.m, d.ϕ
    u = rand(rng, Distributions.InverseGaussian(1, ϕ))
    return rand(rng, Distributions.Poisson(m * u))
end