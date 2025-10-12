"""
    Zipf(N, s)

A *Zipf distribution* ... 

```math
P(X = k) = \\frac{1}{H_{N,s}}\\cdot \\frac{1}{k^s}
```
Where ``H_{N, s}`` is a generalized harmonic number

```julia
Zipf()      # equivalent to BetaNegBinomial(1, 1)
Zipf(N)     # equivalent to BetaNegBinomial(N, 1)

params(d)   # Get the parameters, i.e. (N, s)
```

External link:

*[Zipf distribution on Wikipedia](https://en.wikipedia.org/wiki/Zipf%27s_law)
"""
struct Zipf{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    N::Int
    s::T
    Zipf{T}(N, s) where {T<:Real} = new{T}(N, s)
end

Zipf(N::Integer, s::Integer; check_args::Bool=true) = Zipf(N, float(s); check_args=check_args)
function Zipf(N::Integer, s::Real; check_args::Bool=true)
    @check_args Zipf (N, N > zero(N)) (s, zero(s) <= s)
    return Zipf{typeof(s)}(N, s)
end

Zipf() = Zipf{Float64}(1, 1)
Zipf(N::Integer) = Zipf{Float64}(N, 1)

@distr_support Zipf 1 d.N

#function Distributions.insupport(d::Zipf, x::Int)
#    return 1 <= x <= d.N
#end

# parameters
params(d::Zipf) = (d.N, d.s)
# statistics 
Statistics.mean(d::Zipf) = harmonic(d.N, d.s - 1)/harmonic(d.N, d.s)
Statistics.var(d::Zipf) = harmonic(d.N, d.s - 2)/harmonic(d.N, d.s) - (harmonic(d.N, d.s - 1)/harmonic(d.N, d.s))^2
#Statistics.median(d::Zipf) = 
StatsBase.mode(d::Zipf) = 1 

#evaluate functions CDF, PDF, logPDF, Quantil
function Distributions.cdf(d::Zipf, x::Real)
    N, s = d.N, d.s
    if x < 1
        return 0.0
    elseif x >= N
        return 1.0
    end
    x = round(Int, x)
    num = harmonic(x, s)
    den = harmonic(N, s)
    return num / den
end

function Distributions.pdf(d::Zipf, x::Real)
    N, s = d.N, d.s
    if !insupport(d, x)
        return 0.0
    end
    x = round(Int, x)
    num = 1 / x^s
    den = harmonic(N, s)
    return num / den
end


function Distributions.logpdf(d::Zipf, x::Real)
    N, s = d.N, d.s
    if !insupport(d, x)
        return -Inf
    end
    x = round(Int, x)
    term1 = -s * log(x)
    term2 = log(harmonic(N, s))
    return term1 - term2
end

function Distributions.quantile(d::Zipf, p::Real)
    N, s = d.N, d.s
    if p < 0.0 || p > 1.0
        throw(DomainError(p, "p must be between 0 and 1"))
    end

    if p == 0.0
        return 1
    elseif p == 1.0
        return N
    end

    low, high = 1, N
    while low < high
        mid = div(low + high, 2)
        if cdf(d, mid) < p
            low = mid + 1
        else
            high = mid
        end
    end
    return low
end

# sampling
function Distributions.rand(rng::Distributions.AbstractRNG, d::Zipf)
    u = rand(rng)
    return Distributions.quantile(d, u)
end

function harmonic(N::Integer, s::Real)
    if s == 1
        if N <= 10
            value = 0.0
            for k=1:N
                value +=  1.0 / k
            end
            return value
        else
            return Base.MathConstants.eulergamma + SpecialFunctions.digamma(N+1) # ψ(n+1) + γ = Hₙ 
        end
    else
        if N <= 10 
            value = 0.0
            for k=1:N
                value +=  1.0 / (k^s)
            end
            return value
        else
            return SpecialFunctions.zeta(s) - 1.0 / ((s - 1) * N^(s - 1))
        end
    end
end