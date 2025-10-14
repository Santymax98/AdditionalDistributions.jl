struct MvGaussian{T<:Real, M<:AbstractMatrix{T}, V<:AbstractVector{T}} <: Distributions.ContinuousMultivariateDistribution
    μ::V
    Σ::M
end

MvGaussian(μ::AbstractVector, Σ::AbstractMatrix) = MvGaussian{eltype(μ), typeof(Σ), typeof(μ)}(μ, Σ)
