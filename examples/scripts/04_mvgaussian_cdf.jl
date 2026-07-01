using AdditionalDistributions
using LinearAlgebra
using Random

rng = MersenneTwister(123)

μ = zeros(4)

Σ = [
    1.0 0.6 0.3 0.2
    0.6 1.0 0.5 0.3
    0.3 0.5 1.0 0.4
    0.2 0.3 0.4 1.0
]

d = MvGaussian(μ, Σ)

lower = fill(-1.0, 4)
upper = fill(1.0, 4)

res = cdf_result(
    d,
    lower,
    upper;
    rng=rng,
    m=100_000,
    abseps=1e-5,
    releps=1e-5,
)

println("Distribution:")
println(d)

println("\nRectangular probability:")
println("P(-1 ≤ Xᵢ ≤ 1, i = 1,…,4)")

println("\nCDFResult:")
println("value     = ", res.value)
println("error     = ", res.error)
println("inform    = ", res.inform)
println("neval     = ", res.neval)
println("algorithm = ", res.algorithm)