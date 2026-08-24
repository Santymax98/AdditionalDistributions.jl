"""
    Rademacher()

Rademacher distribution with probability `1/2` at `-1` and probability `1/2` at `1`.
"""
struct Rademacher <: Distributions.DiscreteUnivariateDistribution
end

@distr_support Rademacher -1 1

params(::Rademacher) = ()
Statistics.mean(::Rademacher) = 0
Statistics.var(::Rademacher) = 1
Statistics.median(::Rademacher) = 0
StatsBase.mode(::Rademacher) = NaN
StatsBase.skewness(::Rademacher) = 0
StatsBase.kurtosis(::Rademacher) = -2
StatsBase.entropy(::Rademacher) = log(2)

function Distributions.cdf(::Rademacher, x::Real)
    isnan(float(x)) && return NaN
    x < -1 && return 0.0
    x < 1 && return 0.5
    return 1.0
end

function Distributions.pdf(::Rademacher, x::Real)
    return (x == -1 || x == 1) ? 0.5 : 0.0
end

function Distributions.logpdf(::Rademacher, x::Real)
    return (x == -1 || x == 1) ? log(0.5) : -Inf
end

Distributions.mgf(::Rademacher, t) = cosh(t)
Distributions.cf(::Rademacher, t) = cos(t)

function Distributions.quantile(::Rademacher, p::Real)
    (0 <= p <= 1) || throw(DomainError(p, "p must be in [0, 1]"))
    return p <= 0.5 ? -1 : 1
end

function Distributions.rand(rng::Distributions.AbstractRNG, ::Rademacher)
    return rand(rng, Bool) ? 1 : -1
end
