using AdditionalDistributions
using BenchmarkTools
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using MvNormalCDF
using Random

const SEED = 1234

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

function run_additional(μ, Σ, a, b, maxpts, abseps, releps, seed)
    dist = MvGaussian(μ, Σ)
    rng = MersenneTwister(seed)
    return cdf_result(dist, a, b;
        m=maxpts,
        abseps=abseps,
        releps=releps,
        rng=rng
    )
end

function run_mvnormalcdf(μ, Σ, a, b, maxpts, seed)
    rng = MersenneTwister(seed)
    p, e = MvNormalCDF.mvnormcdf(μ, Σ, a, b; m=maxpts, rng=rng)
    return (value=p, error=e, inform=0)
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

    maxpts = 100_000
    abseps = 1e-6
    releps = 1e-6

    ref = case.reference == :exact ? exact_independent_normal(μ, Σ, a, b) : missing

    # AdditionalDistributions.jl
    res_add = run_additional(μ, Σ, a, b, maxpts, abseps, releps, SEED)

    bench_add = @benchmark run_additional($μ, $Σ, $a, $b, $maxpts, $abseps, $releps, $SEED) samples=5 evals=1

    abs_err_add = ismissing(ref) ? missing : abs(res_add.value - ref)

    push!(rows, (
        "AdditionalDistributions",
        "MvGaussian",
        d,
        case.matrix_case,
        maxpts,
        abseps,
        releps,
        Float64(res_add.value),
        Float64(res_add.error),
        ref,
        abs_err_add,
        median(bench_add.times) / 1e9,
        median(bench_add.allocs),
        median(bench_add.memory),
        res_add.inform,
        String(res_add.algorithm),
        SEED
    ))

    # MvNormalCDF.jl
    res_mvn = run_mvnormalcdf(μ, Σ, a, b, maxpts, SEED)

    bench_mvn = @benchmark run_mvnormalcdf($μ, $Σ, $a, $b, $maxpts, $SEED) samples=5 evals=1

    abs_err_mvn = ismissing(ref) ? missing : abs(res_mvn.value - ref)

    push!(rows, (
        "MvNormalCDF",
        "MvGaussian",
        d,
        case.matrix_case,
        maxpts,
        abseps,
        releps,
        Float64(res_mvn.value),
        Float64(res_mvn.error),
        ref,
        abs_err_mvn,
        median(bench_mvn.times) / 1e9,
        median(bench_mvn.allocs),
        median(bench_mvn.memory),
        res_mvn.inform,
        "qsimvnv",
        SEED
    ))
end

CSV.write(joinpath(@__DIR__, "results", "mvn_compare_mvnormalcdf.csv"), rows)

println(rows)
println()
println("Saved benchmark/results/mvn_compare_mvnormalcdf.csv")