"""
    Conway(λ, ν)

A *Conway–Maxwell–Poisson* distribution, often used to model overdispersed and underdispersed count data, is defined by the following probability mass function (PMF):

```math
P(X = x) = \\frac{\\lambda^x}{(x!)^\\nu Z(\\lambda, \\nu)}, \\quad x \\in \\{0, 1, 2, \\dots\\}
```
where:

- ``Z(\\lambda, \\nu) = \\sum_{j=0}^{\\infty} \\frac{\\lambda^j}{(j!)^\\nu}`` is a normalization constant that ensures the sum of probabilities equals 1.

```julia
Conway()        # equivalent to Conway(1, 1)
Conway(λ)       # equivalent to Conway(λ, 1)
Conway(λ, ν)    # equivalent to Conway(λ, ν)

params(d)        # Get the parameters, i.e. (λ, ν)
```

External link:

* [Conway Maxwell Poisson distribution on Wikipedia](https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution)
"""
struct Conway{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    λ::T
    ν::T
    Conway{T}(λ, ν) where {T<:Real} = new{T}(λ, ν)
end

# Constructor for integer inputs
Conway(λ::Integer, ν::Integer; check_args::Bool=true) = Conway(float(λ), float(ν); check_args=check_args)

# Main constructor with argument checking
function Conway(λ::Real, ν::Real; check_args::Bool=true)
    @check_args Conway (λ, zero(λ) < λ) (ν, zero(ν) <= ν)
    return Conway{typeof(λ)}(λ, ν)
end

Conway() = Conway{Float64}(1, 1)
Conway(λ) = Conway{Float64}(λ, 1.0)

@distr_support Conway 0 Inf
# Parameters
params(d::Conway) = (d.λ, d.ν)

Base.eltype(::Type{Conway{T}}) where {T} = T

# Normalization constant
function Ζ(d::Conway; tol=1e-15, max_iterations=10000)
    λ, ν = d.λ, d.ν
    if ν == 0
        return (1 - λ)^(-1)
    elseif ν == 1
        return exp(λ)
    elseif ν == 2
        return SpecialFunctions.besseli(0, 2 * sqrt(λ))
    else
        Z = 0.0
        previous_term = 1.0
        k = 0
        while k < max_iterations
            log_term = k * log(λ) - ν * SpecialFunctions.loggamma(k + 1)
            term = exp(log_term)
            if term / previous_term < tol
                break
            end
            Z += term
            previous_term = term
            k += 1
        end
        return Z
    end
end

# Some statistics
function Statistics.mean(d::Conway; tol=1e-15, max_iterations=10000)
    λ, ν = d.λ, d.ν
    Ζ_d = Ζ(d; tol=tol, max_iterations=max_iterations)
    log_Z_d = log(Ζ_d)
    μ = 0.0
    previous_term = 1.0
    k = 1  # log(0) non defined
    while k < max_iterations
        log_term = log(k) + k * log(λ) - ν * SpecialFunctions.loggamma(k + 1) - log_Z_d
        term = exp(log_term)
        if term / previous_term < tol
            break
        end
        μ += term
        previous_term = term
        k += 1
    end
    return μ
end

function Statistics.var(d::Conway; tol=1e-15, max_iterations=10000)
    λ, ν = d.λ, d.ν
    Ζ_d = Ζ(d; tol=tol, max_iterations=max_iterations)
    log_Z_d = log(Ζ_d)
    μ = mean(d; tol=tol, max_iterations=max_iterations)
    σ² = 0.0
    previous_term = 1.0
    k = 1  # log(0) non defined
    while k < max_iterations
        log_term = log(k^2) + k * log(λ) - ν * SpecialFunctions.loggamma(k + 1) - log_Z_d
        term = exp(log_term)
        if term / previous_term < tol
            break
        end
        σ² += term
        previous_term = term
        k += 1
    end
    return σ² - μ^2
end

# Evaluate functions CDF, PDF, logPDF, quantile
function Distributions.cdf(d::Conway, x::Real; tol=1e-15, max_iterations=10000)
    λ, ν = d.λ, d.ν
    isnan(float(x)) && return NaN
    x < 0 && return 0.0
    x == Inf && return 1.0
    !isfinite(x) && return 0.0

    x = round(Int, x)
    Z_d = Ζ(d; tol=tol, max_iterations=max_iterations)
    log_Z_d = log(Z_d)
    cdf_value = 0.0
    previous_term = 1.0
    for k in 0:x
        log_term = k * log(λ) - ν * SpecialFunctions.loggamma(k + 1) - log_Z_d
        term = exp(log_term)
        if term / previous_term < tol
            break
        end
        cdf_value += term
        previous_term = term
    end
    return cdf_value
end

function Distributions.pdf(d::Conway, x::Real; tol=1e-15)
    λ, ν = d.λ, d.ν
    if x < 0
        return 0.0
    end
    x = round(Int, x)
    log_Z_d = log(Ζ(d; tol=tol))
    log_term = x * log(λ) - ν * SpecialFunctions.loggamma(x + 1) - log_Z_d
    return exp(log_term)
end

function Distributions.logpdf(d::Conway, x::Real; tol=1e-15)
    λ, ν = d.λ, d.ν
    if x < 0
        return -Inf
    end
    x = round(Int, x)
    log_Z_d = log(Ζ(d; tol=tol))
    return x * log(λ) - ν * SpecialFunctions.loggamma(x + 1) - log_Z_d
end

function Distributions.mgf(d::Conway, t::Real; tol=1e-15)
    λ, ν = d.λ, d.ν
    _d = Conway(exp(t) * λ, ν)
    return Ζ(_d; tol=tol) / Ζ(d; tol=tol)
end

function Distributions.cf(d::Conway, t::Real; tol=1e-15)
    λ, ν = d.λ, d.ν
    _d = Conway(exp(im * t) * λ, ν)
    return Ζ(_d; tol=tol) / Ζ(d; tol=tol)
end

function Distributions.quantile(d::Conway, p::Real; tol=1e-15, max_iterations=10000)
    if p < 0.0 || p > 1.0
        throw(DomainError(p, "p must be in [0, 1]"))
    end
    λ, ν = d.λ, d.ν
    Z_d = Ζ(d; tol=tol, max_iterations=max_iterations)
    log_Z_d = log(Z_d)
    cumulative = 0.0
    k = 0
    previous_term = 1.0
    while cumulative < p && k < max_iterations
        log_term = k * log(λ) - ν * SpecialFunctions.loggamma(k + 1) - log_Z_d
        term = exp(log_term)
        if term / previous_term < tol
            break
        end
        cumulative += term
        previous_term = term
        if cumulative >= p
            return k
        end
        k += 1
    end
    return k
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Conway)
    u = rand(rng)
    return Distributions.quantile(d, u)
end