"""
    ZIP(λ, p)

The *Zero-Inflated Poisson (ZIP) distribution* is a discrete probability distribution that models a scenario where there are more zeros in the data than would be expected from a standard Poisson distribution. It is defined by the following probability mass function (PMF):

```math
P(X = k) =
\\begin{cases} 
p + (1 - p) \\cdot e^{-\\lambda} & \\text{if } k = 0, \\\\
(1 - p) \\cdot \\frac{\\lambda^k \\cdot e^{-\\lambda}}{k!} & \\text{if } k \\geq 1.
\\end{cases}

```julia
ZIP()      # equivalent to ZIP(1, 0.5)
ZIP(λ)     # equivalent to ZIP(λ, 0.5)

params(d)   # Get the parameters, i.e. (λ, p)
```

External link:

*[Zero Inflated Poisson (ZIP) distribution on Wikipedia](https://en.wikipedia.org/wiki/Zero-inflated_model#Zero-inflated_Poisson)
"""
struct ZIP{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    λ::T
    p::T
    ZIP{T}(λ, p) where {T<:Real} = new{T}(λ, p)
end

ZIP(λ::Integer, p::Integer; check_args::Bool=true) = ZIP(float(λ), float(p); check_args=check_args)
function ZIP(λ::Real, p::Real; check_args::Bool=true)
    @check_args ZIP (λ, λ > zero(λ)) (p, zero(p) <= p <= one(p))
    return ZIP{typeof(λ)}(λ, p)
end

ZIP() = ZIP{Float64}(1.0, 0.5)
ZIP(λ) = ZIP{Float64}(λ, 0.5)

@distr_support ZIP 0 Inf


# parameters
params(d::ZIP) = (d.λ, d.p)
# statistics 
Statistics.mean(d::ZIP) = d.λ * (1 - d.p)
Statistics.var(d::ZIP) = d.λ * (1 - d.p) * (1 + d.p * d.λ)
#evaluate functions CDF, PDF, logPDF, Quantil
function Distributions.cdf(d::ZIP, x::Real)
    λ, p = d.λ, d.p
    isnan(float(x)) && return NaN
    x < 0 && return 0.0
    x == Inf && return 1.0
    !isfinite(x) && return 0.0

    x = floor(Int, x)
    cdf_value = 0.0
    for k in 0:x
        if k == 0
            cdf_value += p + (1 - p) * exp(-λ)
        else
            log_pdf_k = LogExpFunctions.log1p(-p) + k * log(λ) - λ - SpecialFunctions.logfactorial(k)
            cdf_value += exp(log_pdf_k)
        end
    end
    return cdf_value
end

# PDF function
function Distributions.pdf(d::ZIP, x::Real)
    λ, p = d.λ, d.p
    if !insupport(d, x)
        return 0.0
    end
    x = round(Int, x)
    if x == 0
        pdf_value = p + (1 - p) * exp(-λ)
    else
        log_pdf_value = log1p(-p) + x * log(λ) - λ - SpecialFunctions.logfactorial(x)
        pdf_value = exp(log_pdf_value)
    end
    return pdf_value
end

# logPDF function
function Distributions.logpdf(d::ZIP, x::Real)
    λ, p = d.λ, d.p
    if !insupport(d, x)
        return -Inf
    end
    x = round(Int, x)
    if x == 0
        log_pdf_value = log(p + (1 - p) * exp(-λ))
    else
        log_pdf_value = log1p(-p) + x * log(λ) - λ - SpecialFunctions.logfactorial(x)
    end
    return log_pdf_value
end

function Distributions.quantile(d::ZIP, p::Real)
    @assert 0.0 <= p <= 1.0 "p must be between 0 and 1"
    
    λ, pi = d.λ, d.p
    
    # Define the initial search range for the quantile
    lo, hi = 0, max(1, floor(Int, λ))
    
    # Expand the range until the CDF at `hi` is greater than or equal to `p`
    while Distributions.cdf(d, hi) < p
        lo = hi
        hi *= 2
    end
    
    # Binary search for the quantile
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

# sampling
function Distributions.rand(rng::Distributions.AbstractRNG, d::ZIP)
    λ, p = d.λ, d.p
    u = rand(rng)
    if u < p
        return 0
    else
        return rand(rng, Distributions.Poisson(λ))
    end
end