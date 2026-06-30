@testitem "Multivariate CDF invariants – MvGaussian" tags=[:multivariate] begin
    using LinearAlgebra
    using Distributions
    using Statistics
    using AdditionalDistributions

    @testset "d=1 agrees with Normal" begin
        μ = [1.0]
        Σ = reshape([4.0], 1, 1)
        d = MvGaussian(μ, Σ)

        a = [-0.5]
        b = [2.0]

        ref = cdf(Normal(1.0, 2.0), b[1]) - cdf(Normal(1.0, 2.0), a[1])
        val = cdf(d, a, b; m=50_000, abseps=1e-8, releps=1e-8)

        @test val ≈ ref atol=1e-8
    end

    @testset "diagonal covariance factorizes" begin
        μ = [0.0, 1.0, -1.0]
        σ = [1.0, 2.0, 0.5]
        Σ = Matrix(Diagonal(σ .^ 2))

        d = MvGaussian(μ, Σ)

        a = [-1.0, -2.0, -1.5]
        b = [1.0, 3.0, 0.5]

        ref = prod(
            cdf(Normal(μ[i], σ[i]), b[i]) -
            cdf(Normal(μ[i], σ[i]), a[i])
            for i in eachindex(μ)
        )

        val = cdf(d, a, b; m=100_000, abseps=1e-7, releps=1e-7)

        @test val ≈ ref atol=5e-5
    end

    @testset "total probability" begin
        μ = zeros(3)
        Σ = Matrix(I, 3, 3)
        d = MvGaussian(μ, Σ)

        a = fill(-Inf, 3)
        b = fill(Inf, 3)

        @test cdf(d, a, b) ≈ 1.0 atol=1e-12
    end

    @testset "empty rectangle" begin
        μ = zeros(2)
        Σ = Matrix(I, 2, 2)
        d = MvGaussian(μ, Σ)

        a = [1.0, 0.0]
        b = [0.0, 1.0]

        @test cdf(d, a, b) == 0.0
    end
end


@testitem "Multivariate CDF invariants – MvTStudent" tags=[:multivariate] begin
    using LinearAlgebra
    using Distributions
    using AdditionalDistributions

    @testset "empty rectangle" begin
        μ = zeros(2)
        Σ = Matrix(I, 2, 2)
        d = MvTStudent(4.5, μ, Σ)

        a = [1.0, 0.0]
        b = [0.0, 1.0]

        @test cdf(d, a, b) == 0.0
    end

    @testset "non-integer degrees of freedom" begin
        μ = zeros(2)
        Σ = Matrix(I, 2, 2)
        d = MvTStudent(2.5, μ, Σ)

        a = [-1.0, -1.0]
        b = [1.0, 1.0]

        val = cdf(d, a, b; m=80_000, abseps=1e-6, releps=1e-6)

        @test isfinite(val)
        @test 0.0 <= val <= 1.0
    end
end