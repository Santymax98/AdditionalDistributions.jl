using AdditionalDistributions
using Distributions
using Random
using Statistics

rng = MersenneTwister(123)

examples = [
    ("Lomax", Lomax(2.5, 1.0), 1.0),
    ("Burr", Burr(2.0, 3.0), 1.0),
    ("Dagum", Dagum(2.0, 3.0), 1.0),
    ("Maxwell", Maxwell(1.5), 1.0),
    ("Gompertz", Gompertz(1.2, 0.8), 1.0),
    ("PERT", PERT(0.0, 0.6, 1.0), 0.5),
]

println("Basic Distributions.jl interface examples")
println("=========================================")

for (name, d, x0) in examples
    sample = rand(rng, d, 5_000)

    println("\n$name")
    println("-"^length(name))
    println("distribution       = ", d)
    println("pdf(d, x0)         = ", pdf(d, x0))
    println("logpdf(d, x0)      = ", logpdf(d, x0))
    println("cdf(d, x0)         = ", cdf(d, x0))
    println("quantile(d, 0.95)  = ", quantile(d, 0.95))
    println("empirical mean     = ", mean(sample))

    try
        println("theoretical mean   = ", mean(d))
    catch err
        println("theoretical mean   = not available")
    end
end