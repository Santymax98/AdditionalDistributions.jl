using AdditionalDistributions
using LinearAlgebra
using Random

rng = MersenneTwister(123)

ν = 4.0
μ = zeros(3)

Σ = [
    1.0 0.7 0.5
    0.7 1.0 0.6
    0.5 0.6 1.0
]

d = MvTStudent(ν, μ, Σ)

# Probability of a joint lower-tail event.
# This can be interpreted as the probability that three dependent
# heavy-tailed variables are simultaneously between -2 and 0.
lower = fill(-2.0, 3)
upper = fill(0.0, 3)

res = cdf_result(
    d,
    lower,
    upper;
    rng=rng,
    m=200_000,
    abseps=2e-5,
    releps=2e-5,
)

println("Distribution:")
println(d)

println("\nRectangular probability:")
println("P(-2 ≤ Xᵢ ≤ 0, i = 1,…,3)")

println("\nCDFResult:")
println("value     = ", res.value)
println("error     = ", res.error)
println("inform    = ", res.inform)
println("neval     = ", res.neval)
println("algorithm = ", res.algorithm)