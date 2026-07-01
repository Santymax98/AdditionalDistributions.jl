using AdditionalDistributions
using Distributions
using Random
using Statistics

rng = MersenneTwister(123)

models = [
    ("Lomax", Lomax(2.2, 1_000.0)),
    ("Burr", Burr(1.6, 2.5, 1_000.0)),
    ("Dagum", Dagum(2.0, 1_000.0, 1.5)),
]

println("Heavy-tailed loss models")
println("========================")

for (name, d) in models
    losses = rand(rng, d, 100_000)

    VaR_95 = quantile(losses, 0.95)
    VaR_99 = quantile(losses, 0.99)

    ES_95 = mean(losses[losses .>= VaR_95])
    ES_99 = mean(losses[losses .>= VaR_99])

    println("\n$name")
    println("-"^length(name))
    println("distribution = ", d)
    println("mean         = ", mean(losses))
    println("median       = ", median(losses))
    println("VaR 95%      = ", VaR_95)
    println("VaR 99%      = ", VaR_99)
    println("ES 95%       = ", ES_95)
    println("ES 99%       = ", ES_99)
    println("theoretical q95 = ", quantile(d, 0.95))
    println("theoretical q99 = ", quantile(d, 0.99))
end