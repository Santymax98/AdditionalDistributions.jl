using TestItems

@testitem "MvTStudent – native wrapper interface" tags=[:multivariate, :api, :student] begin
    using Test
    using AdditionalDistributions
    using Distributions
    using Statistics
    using Random

    ν = 5.0
    μ = [0.5, -1.0]
    Σ = [1.0 0.25;
         0.25 2.0]

    d  = MvTStudent(ν, μ, Σ)
    dn = MvTDist(ν, μ, Σ)

    x = [0.2, -0.3]

    @test d.dist isa Distributions.AbstractMvTDist
    @test length(d) == length(dn)

    pd = params(d)
    pn = params(dn)
    @test pd[1] == pn[1]
    @test pd[2] ≈ pn[2]
    @test Matrix(pd[3]) ≈ Matrix(pn[3])

    @test mean(d) ≈ mean(dn)
    @test var(d) ≈ var(dn)
    @test cov(d) ≈ cov(dn)
    @test scale(d) ≈ scale(dn)

    @test pdf(d, x) ≈ pdf(dn, x)
    @test logpdf(d, x) ≈ logpdf(dn, x)
    @test entropy(d) ≈ entropy(dn)
    @test insupport(d, x) == insupport(dn, x)

    rng1 = MersenneTwister(2024)
    rng2 = MersenneTwister(2024)
    @test rand(rng1, d) ≈ rand(rng2, dn)
end


@testitem "MvTStudent – native CDF backend" tags=[:multivariate, :api, :student, :cdf] begin
    using Test
    using AdditionalDistributions
    using Distributions

    ν = 5.0
    μ = [2.0]
    Σ = reshape([9.0], 1, 1)

    d  = MvTStudent(ν, μ, Σ)
    dn = MvTDist(ν, μ, Σ)

    a = [-Inf]
    b = [2.0]

    # Symmetry about the location gives an exact 1/2 probability,
    # independently of ν and the scale.
    expected = 0.5

    rn = cdf_result(dn, a, b)
    rw = cdf_result(d, a, b)

    @test rn.value ≈ expected atol=1e-14
    @test rw.value ≈ expected atol=1e-14
    @test rn.algorithm == :univariate_exact_t
    @test rw.algorithm == :univariate_exact_t

    @test AdditionalDistributions._cdf(dn, a, b) ≈ expected atol=1e-14
    @test AdditionalDistributions._cdf(dn, b) ≈ expected atol=1e-14

    @test cdf(d, a, b) ≈ expected atol=1e-14
    @test cdf(d, b) ≈ expected atol=1e-14

    @test_throws DimensionMismatch cdf(d, [0.0, 1.0])
    @test_throws DimensionMismatch cdf_result(dn, [-Inf, -Inf], b)
end


@testitem "MvTStudent – CDF uses location and scale" tags=[:multivariate, :api, :student, :cdf] begin
    using Test
    using AdditionalDistributions
    using Distributions

    # ν < 1: the statistical mean is undefined, but the CDF is perfectly
    # well-defined around the location μ.
    ν = 0.5
    μ = [3.0]
    Σ = reshape([4.0], 1, 1)

    d  = MvTStudent(ν, μ, Σ)
    dn = MvTDist(ν, μ, Σ)

    @test !(mean(dn) isa AbstractVector)

    @test cdf_result(dn, [-Inf], μ).value ≈ 0.5 atol=1e-14
    @test cdf(d, μ) ≈ 0.5 atol=1e-14
end
