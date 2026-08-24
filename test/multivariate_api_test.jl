using TestItems

@testitem "MvGaussian – AbstractMvNormal interface" tags=[:multivariate, :api] begin
    using Test
    using AdditionalDistributions
    using Distributions
    using Statistics
    using Random
    using LinearAlgebra

    μ = [0.0, 1.0]
    Σ = [1.0 0.4;
         0.4 2.0]

    d  = MvGaussian(μ, Σ)
    dn = MvNormal(μ, Σ)

    x = [0.3, -0.2]

    @test d isa Distributions.AbstractMvNormal

    @test length(d) == length(dn)
    @test mean(d) == mean(dn)
    @test var(d) ≈ var(dn)
    @test Matrix(cov(d)) ≈ Matrix(cov(dn))

    @test Matrix(invcov(d)) ≈ Matrix(invcov(dn))
    @test logdetcov(d) ≈ logdetcov(dn)
    @test sqmahal(d, x) ≈ sqmahal(dn, x)
    @test gradlogpdf(d, x) ≈ gradlogpdf(dn, x)

    @test pdf(d, x) ≈ pdf(dn, x)
    @test logpdf(d, x) ≈ logpdf(dn, x)
    @test entropy(d) ≈ entropy(dn)

    rng1 = MersenneTwister(2024)
    rng2 = MersenneTwister(2024)
    @test rand(rng1, d) ≈ rand(rng2, dn)

    @test insupport(d, x) == insupport(dn, x)
end


@testitem "MvGaussian – native CDF backend" tags=[:multivariate, :api, :cdf] begin
    using Test
    using AdditionalDistributions
    using Distributions

    μ = [0.0, 1.0]
    Σ = [1.0 0.0;
         0.0 4.0]

    d  = MvGaussian(μ, Σ)
    dn = MvNormal(μ, Σ)

    a = [-Inf, -Inf]
    b = [0.0, 1.0]

    expected = 0.25

    rn = cdf_result(dn, a, b)
    rw = cdf_result(d, a, b)

    @test rn.value ≈ expected atol=1e-14
    @test rw.value ≈ expected atol=1e-14
    @test rn.value ≈ rw.value

    @test AdditionalDistributions._cdf(dn, a, b) ≈ expected atol=1e-14
    @test AdditionalDistributions._cdf(dn, b) ≈ expected atol=1e-14

    @test cdf(d, a, b) ≈ expected atol=1e-14
    @test cdf(d, b) ≈ expected atol=1e-14

    # Native MvNormal receives CDF functionality only through functions
    # owned by AdditionalDistributions; Distributions.cdf is not pirated.
    @test cdf_result(dn, a, b).value == AdditionalDistributions._cdf(dn, a, b)

    @test_throws DimensionMismatch cdf(d, [0.0])
    @test_throws DimensionMismatch cdf_result(dn, [-Inf], b)
end
