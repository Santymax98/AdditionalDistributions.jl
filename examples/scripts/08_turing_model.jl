using AdditionalDistributions
using Distributions
using Random
using Statistics
using Turing

rng = MersenneTwister(123)

println("Turing.jl + AdditionalDistributions.jl")
println("=====================================")

# Simulated positive data from a Maxwell distribution.
true_σ = 1.5
d_true = Maxwell(true_σ)

x = rand(rng, d_true, 200)

@model function maxwell_model(x)
    σ ~ truncated(Normal(1.5, 0.5), 0.1, 5.0)

    for i in eachindex(x)
        x[i] ~ Maxwell(σ)
    end
end

model = maxwell_model(x)

println("\nData summary:")
println("sample size = ", length(x))
println("mean(x)     = ", mean(x))
println("std(x)      = ", std(x))

println("\nModel:")
println(model)

println("\nRunning a short NUTS chain...")
chain = sample(rng, model, NUTS(0.65), 500)

println("\nPosterior summary:")
println(chain)

σ_samples = vec(Array(chain[:σ]))

println("\nPosterior estimate:")
println("true σ          = ", true_σ)
println("posterior mean  = ", mean(σ_samples))
println("posterior std   = ", std(σ_samples))
println("posterior q025  = ", quantile(σ_samples, 0.025))
println("posterior q975  = ", quantile(σ_samples, 0.975))