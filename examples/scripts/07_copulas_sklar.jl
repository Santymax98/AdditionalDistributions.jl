using AdditionalDistributions
using Copulas
using Distributions
using Random
using Statistics
using StatsBase

rng = MersenneTwister(123)

println("Copulas.jl + AdditionalDistributions.jl")
println("=======================================")

# Dependence model: Gaussian copula with positive dependence.
Σ = [
    1.0 0.65
    0.65 1.0
]

C = GaussianCopula(Σ)

# Marginal models from AdditionalDistributions.jl.
# Both are positive heavy-tailed distributions, but with different tail behavior.
marginals = (
    Burr(1.6, 2.5, 1_000.0),
    Dagum(2.0, 1_000.0, 1.5),
)

joint = SklarDist(C, marginals)

n = 5_000
X = rand(rng, joint, n)

println("\nJoint model:")
println(joint)

println("\nSample size:")
println(size(X))

x1 = vec(X[1, :])
x2 = vec(X[2, :])

println("\nMarginal summaries:")
println("mean variable 1   = ", mean(x1))
println("mean variable 2   = ", mean(x2))
println("median variable 1 = ", median(x1))
println("median variable 2 = ", median(x2))
println("q95 variable 1    = ", quantile(x1, 0.95))
println("q95 variable 2    = ", quantile(x2, 0.95))

println("\nDependence summary:")
println("Pearson correlation  = ", cor(x1, x2))
println("Spearman correlation = ", corspearman(x1, x2))

println("\nJoint tail event:")
threshold1 = quantile(marginals[1], 0.95)
threshold2 = quantile(marginals[2], 0.95)

joint_tail_prob = mean((x1 .> threshold1) .& (x2 .> threshold2))

println("threshold 1 = ", threshold1)
println("threshold 2 = ", threshold2)
println("P(X₁ > q₀.₉₅, X₂ > q₀.₉₅) ≈ ", joint_tail_prob)