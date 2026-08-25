@testitem "Fast CBC lattice construction and cache" begin
    using AdditionalDistributions

    const AD = AdditionalDistributions

    α, p = AD._qmc_lattice(Float64, 3, 625)

    @test p == 619
    @test length(α) == 3
    @test round.(Int, α .* p) == [1, 259, 117]
    @test all(x -> 0.0 < x < 1.0, α)

    α2, p2 = AD._qmc_lattice(Float64, 3, 625)

    @test p2 == p
    @test α2 === α

    αsmall, nsmall =
        AD._qmc_lattice(Float64, 3, 2)

    @test nsmall == 2
    @test αsmall ==
        AD.richtmyer_roots(Float64, 4)
end

@testitem "CBC-backed MVN preserves public CDF contract" begin
    using AdditionalDistributions
    using Distributions
    using StableRNGs

    Σ = [
        1.0 0.3 0.3
        0.3 1.0 0.3
        0.3 0.3 1.0
    ]

    d = MvNormal(zeros(3), Σ)
    a = fill(-1.0, 3)
    b = fill(1.0, 3)

    result = AdditionalDistributions.cdf_result(
        d,
        a,
        b;
        m=10_000,
        abseps=0.0,
        releps=0.0,
        nshifts=12,
        rng=StableRNG(20260824),
    )

    @test 0.0 <= result.value <= 1.0
    @test isapprox(
        result.value,
        0.3383461438763211;
        atol=2e-5,
        rtol=0,
    )

    @test result.neval == 10_000
end

@testitem "CBC plus radial tent MVT smoke test" begin
    using AdditionalDistributions
    using StableRNGs

    Σ = [
        1.0 0.3 0.3
        0.3 1.0 0.3
        0.3 0.3 1.0
    ]

    d = MvTStudent(4.0, Σ)
    a = fill(-1.0, 3)
    b = fill(1.0, 3)

    result = cdf_result(
        d,
        a,
        b;
        m=10_000,
        abseps=0.0,
        releps=0.0,
        nshifts=12,
        rng=StableRNG(20260824),
    )

    @test 0.0 <= result.value <= 1.0
    @test isapprox(
        result.value,
        0.31181962814543895;
        atol=5e-4,
        rtol=0,
    )
    @test result.neval == 10_000
end

@testitem "Gaussian default QMC shift count" begin
    using AdditionalDistributions
    using Distributions
    using LinearAlgebra
    using Random

    d = 3
    ρ = 0.3

    Σ = fill(ρ, d, d)
    Σ[diagind(Σ)] .= 1.0

    dist = MvNormal(zeros(d), Σ)
    lower = fill(-1.0, d)
    upper = fill(1.0, d)

    seed = 20260824

    default_result = cdf_result(
        dist,
        lower,
        upper;
        m=10_000,
        rng=MersenneTwister(seed),
    )

    explicit_result = cdf_result(
        dist,
        lower,
        upper;
        m=10_000,
        nshifts=10,
        rng=MersenneTwister(seed),
    )

    @test default_result.value == explicit_result.value
    @test default_result.error == explicit_result.error
    @test default_result.inform == explicit_result.inform
    @test default_result.neval == explicit_result.neval
end
