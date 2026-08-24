"""
    ZIB(n, θ, p)

The *Zero-Inflated Binomial (ZIB) distribution* is a discrete probability distribution that extends the binomial distribution by incorporating an excess of zero counts. The probability mass function (PMF) is defined as:

```math
P(X = k) =
\\begin{cases} 
\\theta + (1 - \\theta) \\cdot (1 - p)^n & \\text{if } k = 0, \\\\
(1 - \\theta) \\cdot \\binom{n}{k} p^k (1 - p)^{n-k} & \\text{if } k > 0.
\\end{cases}
```

```julia
ZIB()       # equivalent to ZIB(1, 0.5, 0.5)
ZIB(n)      # equivalent to ZIB(n, 0.5, 0.5)
ZIB(n, θ)   # equivalent to ZIB(n, θ, 0.5)

params(d)   # Get the parameters, i.e. (n, θ, p)
```

"""
struct ZIB{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    n::Int
    θ::T
    p::T
    ZIB{T}(n, θ, p) where {T<:Real} = new{T}(n, θ, p)
end

ZIB(n::Integer, θ::Integer, p::Integer; check_args::Bool=true) = ZIB(n, float(θ), float(p); check_args=check_args)
function ZIB(n::Integer, θ::Real, p::Real; check_args::Bool=true)
    @check_args ZIB (n, n >= zero(n)) (θ, zero(θ) <= θ <= one(θ)) (p, zero(p) <= p <= one(p))
    return ZIB{typeof(θ)}(n, θ, p)
end

ZIB() = ZIB{Float64}(1, 0.5, 0.5)
ZIB(n) = ZIB{Float64}(n, 0.5, 0.5)
ZIB(n, θ) = ZIB{Float64}(n, θ, 0.5)

@distr_support ZIB 0 d.n

# parameters
params(d::ZIB) = (d.n, d.θ, d.p)

# statistics
Statistics.mean(d::ZIB) = (1 - d.p) * d.n * d.θ
Statistics.var(d::ZIB) = (1 - d.p) * (d.n * d.θ * (1 - d.θ) + d.p * d.n^2 * d.θ^2)

# CDF function
function Distributions.cdf(d::ZIB, x::Real)
    n, θ, p = d.n, d.θ, d.p
    isnan(float(x)) && return NaN
    x < 0 && return 0.0
    x >= n && return 1.0

    x = floor(Int, x)
    cdf_value = 0.0
    for k in 0:x
        if k == 0
            cdf_value += p + (1 - p) * (1 - θ)^n
        else
            log_pdf_k = LogExpFunctions.log1p(-p) + k * log(θ) + (n - k) * LogExpFunctions.log1p(-θ) + SpecialFunctions.logfactorial(n) - SpecialFunctions.logfactorial(k) - SpecialFunctions.logfactorial(n-k)
            cdf_value += exp(log_pdf_k)
        end
    end
    return clamp(cdf_value, 0.0, 1.0)
end

# PDF function
function Distributions.pdf(d::ZIB, x::Real)
    n, θ, p = d.n, d.θ, d.p
    if !insupport(d, x)
        return 0.0
    end
    x = round(Int, x)
    if x == 0
        pdf_value = p + (1 - p) * (1 - θ)^n
    else
        log_pdf_x = LogExpFunctions.log1p(-p) + x * log(θ) + (n - x) * LogExpFunctions.log1p(-θ) + SpecialFunctions.logfactorial(n) - SpecialFunctions.logfactorial(x) - SpecialFunctions.logfactorial(n-x)
        pdf_value = exp(log_pdf_x)
    end
    return pdf_value
end

# logPDF function
function Distributions.logpdf(d::ZIB, x::Real)
    n, θ, p = d.n, d.θ, d.p
    if !insupport(d, x)
        return -Inf
    end
    x = round(Int, x)
    if x == 0
        log_pdf_x = log(p + (1 - p) * (1 - θ)^n)
    else
        log_pdf_x = LogExpFunctions.log1p(-p) + x * log(θ) + (n - x) * LogExpFunctions.log1p(-θ) + SpecialFunctions.logfactorial(n) - SpecialFunctions.logfactorial(x) - SpecialFunctions.logfactorial(n-x)
    end
    return log_pdf_x
end

# Quantile function using binary search
function Distributions.quantile(d::ZIB, p::Real)
    @assert 0.0 <= p <= 1.0 "p must be between 0 and 1"
    
    n, θ, p_zero = d.n, d.θ, d.p
    
    lo, hi = 0, n
    while lo < hi
        mid = (lo + hi) ÷ 2
        if Distributions.cdf(d, mid) < p
            lo = mid + 1
        else
            hi = mid
        end
    end
    
    return lo
end

# Sampling function
function Distributions.rand(rng::Distributions.AbstractRNG, d::ZIB)
    n, θ, p = d.n, d.θ, d.p
    u = rand(rng)
    if u < p
        return 0
    else
        return rand(rng, Distributions.Binomial(n, θ))
    end
end