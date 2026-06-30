using AdditionalDistributions
using BenchmarkTools
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Random

Random.seed!(1234)

function equicorrelation(d, ρ)
    Σ = fill(ρ, d, d)
    Σ[diagind(Σ)] .= 1.0
    return Σ
end

function ar1correlation(d, ρ)
    return [ρ^abs(i-j) for i in 1:d, j in 1:d]
end

function exact_independent_normal(μ, Σ, a, b)
    σ = sqrt.(diag(Σ))
    return prod(
        cdf(Normal(μ[i], σ[i]), b[i]) -
        cdf(Normal(μ[i], σ[i]), a[i])
        for i in eachindex(μ)
    )
end

cases = []

for d in (2, 5, 10, 20)
    push!(cases, (
        dimension=d,
        matrix_case="diagonal",
        μ=zeros(d),
        Σ=Matrix(I, d, d),
        a=fill(-1.0, d),
        b=fill(1.0, d),
        reference=:exact
    ))

    for ρ in (0.2, 0.5, 0.9)
        push!(cases, (
            dimension=d,
            matrix_case="equicorrelated_ρ=$(ρ)",
            μ=zeros(d),
            Σ=equicorrelation(d, ρ),
            a=fill(-1.0, d),
            b=fill(1.0, d),
            reference=:none
        ))

        push!(cases, (
            dimension=d,
            matrix_case="ar1_ρ=$(ρ)",
            μ=zeros(d),
            Σ=ar1correlation(d, ρ),
            a=fill(-1.0, d),
            b=fill(1.0, d),
            reference=:none
        ))
    end
end

rows = DataFrame(
    method=String[],
    distribution=String[],
    dimension=Int[],
    matrix_case=String[],
    maxpts=Int[],
    abseps=Float64[],
    releps=Float64[],
    probability=Float64[],
    reported_error=Float64[],
    reference_probability=Union{Missing, Float64}[],
    absolute_error=Union{Missing, Float64}[],
    time_seconds=Float64[],
    allocations=Int[],
    memory_bytes=Int[],
    inform=Int[],
    algorithm=String[],
    seed=Int[]
)

for case in cases
    d = case.dimension
    μ = case.μ
    Σ = case.Σ
    a = case.a
    b = case.b

    dist = MvGaussian(μ, Σ)

    maxpts = 100_000
    abseps = 1e-6
    releps = 1e-6

    res = cdf_result(dist, a, b; m=maxpts, abseps=abseps, releps=releps)

    bench = @benchmark cdf_result($dist, $a, $b; m=$maxpts, abseps=$abseps, releps=$releps) samples=5 evals=1

    ref = case.reference == :exact ? exact_independent_normal(μ, Σ, a, b) : missing
    abs_err = ismissing(ref) ? missing : abs(res.value - ref)

    push!(rows, (
        "AdditionalDistributions",
        "MvGaussian",
        d,
        case.matrix_case,
        maxpts,
        abseps,
        releps,
        Float64(res.value),
        Float64(res.error),
        ref,
        abs_err,
        median(bench.times) / 1e9,
        median(bench.allocs),
        median(bench.memory),
        res.inform,
        String(res.algorithm),
        1234
    ))
end

CSV.write(joinpath(@__DIR__, "results", "mvn_basic.csv"), rows)

println(rows)
println()
println("Saved benchmark/results/mvn_basic.csv")