using AdditionalDistributions
using Distributions
using Random
using Statistics
using StatsBase

rng = MersenneTwister(123)

models = [
    ("ZIP", ZIP(2.0, 0.35)),
    ("ZINB", ZINB(5, 0.30, 0.45)),
    ("BetaNegBinomial", BetaNegBinomial(10, 2.0, 5.0)),
    ("Delaporte", Delaporte(2.0, 5.0, 0.45)),
    ("Yule", Yule(2.5)),
    ("Zeta", Zeta(2.2)),
]

println("Discrete count models")
println("=====================")

for (name, d) in models
    x = rand(rng, d, 10_000)

    zero_fraction = mean(x .== 0)

    println("\n$name")
    println("-"^length(name))
    println("distribution       = ", d)
    println("minimum            = ", minimum(x))
    println("maximum            = ", maximum(x))
    println("empirical mean     = ", mean(x))
    println("empirical var      = ", var(x))
    println("zero fraction      = ", zero_fraction)

    try
        println("theoretical mean   = ", mean(d))
    catch
        println("theoretical mean   = not available")
    end

    try
        println("theoretical var    = ", var(d))
    catch
        println("theoretical var    = not available")
    end

    println("first counts       = ", countmap(x[1:50]))
end