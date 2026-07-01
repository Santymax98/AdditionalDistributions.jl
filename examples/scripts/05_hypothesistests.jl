using AdditionalDistributions
using HypothesisTests
using Distributions
using LinearAlgebra
using Random
using Statistics

rng = MersenneTwister(123)

println("HypothesisTests.jl workflows")
println("============================")

# -------------------------------------------------------------------
# 1. Univariate two-sample comparison
# -------------------------------------------------------------------

println("\n1. Univariate heavy-tailed comparison")
println("------------------------------------")

x = rand(rng, Burr(1.6, 2.5, 1_000.0), 500)
y = rand(rng, Dagum(2.0, 1_000.0, 1.5), 500)

ks = ApproximateTwoSampleKSTest(x, y)
mw = MannWhitneyUTest(x, y)

println("Sample summaries:")
println("mean(x)   = ", mean(x))
println("mean(y)   = ", mean(y))
println("median(x) = ", median(x))
println("median(y) = ", median(y))
println("q95(x)    = ", quantile(x, 0.95))
println("q95(y)    = ", quantile(y, 0.95))

println("\nTwo-sample Kolmogorov-Smirnov test:")
println(ks)

println("\nMann-Whitney U test:")
println(mw)

# -------------------------------------------------------------------
# 2. Multivariate tests
# -------------------------------------------------------------------

println("\n\n2. Multivariate tests")
println("--------------------")

μ1 = [0.0, 0.0, 0.0]
μ2 = [0.35, 0.00, 0.25]

Σ1 = [
    1.0 0.5 0.3
    0.5 1.0 0.4
    0.3 0.4 1.0
]

Σ2 = [
    1.2 0.4 0.2
    0.4 1.1 0.3
    0.2 0.3 1.0
]

d1 = MvGaussian(μ1, Σ1)
d2 = MvGaussian(μ2, Σ2)

n1 = 250
n2 = 250

# Temporary workaround until MvGaussian exposes rand(d, n) directly.
# HypothesisTests.jl expects observations in rows, so we transpose.
X = permutedims(rand(rng, d1, n1))
Y = permutedims(rand(rng, d2, n2))

println("Data dimensions:")
println("size(X) = ", size(X))
println("size(Y) = ", size(Y))

println("\nSample mean vectors:")
println("mean(X, dims=1) = ", vec(mean(X; dims=1)))
println("mean(Y, dims=1) = ", vec(mean(Y; dims=1)))

println("\nUnequal covariance Hotelling T² test:")
println(UnequalCovHotellingT2Test(X, Y))

println("\nEqual covariance Hotelling T² test:")
println(EqualCovHotellingT2Test(X, Y))

println("\nBartlett test for equality of covariance matrices:")
println(BartlettTest(X, Y))