@testitem "Reference values – continuous univariate distributions" tags=[:continuous, :reference] begin
    using AdditionalDistributions
    using Distributions
    using LogExpFunctions
    using SpecialFunctions
    using Statistics
    using Test

    @testset "Burr" begin
        c, k, λ = 2.0, 3.0, 4.0
        d = Burr(c, k, λ)
        x = 2.0
        p = 0.25

        cdf_ref = 1 - (1 + (x / λ)^c)^(-k)
        pdf_ref = (c * k / λ) * (x / λ)^(c - 1) * (1 + (x / λ)^c)^(-(k + 1))
        q_ref = λ * ((1 / (1 - p)^(1 / k)) - 1)^(1 / c)

        @test params(d) == (c, k, λ)
        @test cdf(d, x) ≈ cdf_ref atol=1e-14
        @test pdf(d, x) ≈ pdf_ref atol=1e-14
        @test logpdf(d, x) ≈ log(pdf_ref) atol=1e-14
        @test quantile(d, p) ≈ q_ref atol=1e-14
        @test cdf(d, -1.0) == 0.0
    end

    @testset "Bradford" begin
        a = 5.0
        d = Bradford(a)
        x = 0.4
        p = 0.7

        cdf_ref = log1p(a * x) / log1p(a)
        pdf_ref = a / (log1p(a) * (1 + a * x))
        q_ref = expm1(p * log1p(a)) / a

        @test params(d) == (a,)
        @test cdf(d, x) ≈ cdf_ref atol=1e-14
        @test pdf(d, x) ≈ pdf_ref atol=1e-14
        @test logpdf(d, x) ≈ log(pdf_ref) atol=1e-14
        @test quantile(d, p) ≈ q_ref atol=1e-14
        @test cdf(d, -0.1) == 0.0
        @test cdf(d, 1.0) ≈ 1.0 atol=1e-14
    end

    @testset "Lomax" begin
        α, λ = 2.5, 3.0
        d = Lomax(α, λ)
        x = 1.25
        p = 0.6

        cdf_ref = 1 - (1 + x / λ)^(-α)
        pdf_ref = (α / λ) * (1 + x / λ)^(-(α + 1))
        q_ref = λ * ((1 - p)^(-1 / α) - 1)

        @test params(d) == (α, λ)
        @test cdf(d, x) ≈ cdf_ref atol=1e-14
        @test pdf(d, x) ≈ pdf_ref atol=1e-14
        @test logpdf(d, x) ≈ log(pdf_ref) atol=1e-14
        @test quantile(d, p) ≈ q_ref atol=1e-14
        @test quantile(d, 0.0) == 0.0
        @test quantile(d, 1.0) == Inf
    end

    @testset "Maxwell" begin
        a = 2.0
        d = Maxwell(a)
        x = 1.5

        cdf_ref = erf(x / (sqrt(2) * a)) - sqrt(2 / π) * (x / a) * exp(-(x^2 / (2 * a^2)))
        pdf_ref = sqrt(2 / π) * (x^2 / a^3) * exp(-(x^2 / (2 * a^2)))

        @test params(d) == (a,)
        @test cdf(d, x) ≈ cdf_ref atol=1e-14
        @test pdf(d, x) ≈ pdf_ref atol=1e-14
        @test logpdf(d, x) ≈ log(pdf_ref) atol=1e-14
        @test mean(d) ≈ 2a * sqrt(2 / π) atol=1e-14
        @test cdf(d, -1.0) == 0.0
    end

    @testset "Nakagami" begin
        m, Ω = 2.0, 3.0
        d = Nakagami(m, Ω)
        x = 1.2

        cdf_ref = gamma_inc(m, (m / Ω) * x^2)[1]
        pdf_ref = (2 * m^m / (gamma(m) * Ω^m)) * x^(2m - 1) * exp(-(m / Ω) * x^2)

        @test params(d) == (m, Ω)
        @test cdf(d, x) ≈ cdf_ref atol=1e-14
        @test pdf(d, x) ≈ pdf_ref atol=1e-14
        @test logpdf(d, x) ≈ log(pdf_ref) atol=1e-14
        @test cdf(d, -0.5) == 0.0
    end

    @testset "PERT" begin
        a, b, m = 0.0, 1.0, 3.0
        d = PERT(a, b, m)
        x = 1.2
        z = (x - a) / (m - a)
        α = 1 + 4 * (b - a) / (m - a)
        β = 1 + 4 * (m - b) / (m - a)

        cdf_ref = beta_inc(α, β, z)[1]
        pdf_ref = ((x - a)^(α - 1) * (m - x)^(β - 1)) / (beta(α, β) * (m - a)^(α + β - 1))

        @test params(d) == (a, b, m)
        @test cdf(d, x) ≈ cdf_ref atol=1e-14
        @test pdf(d, x) ≈ pdf_ref atol=1e-14
        @test logpdf(d, x) ≈ log(pdf_ref) atol=1e-14
        @test quantile(d, 0.0) == a
        @test quantile(d, 1.0) == m
    end
end

@testitem "Reference values – discrete univariate distributions" tags=[:discrete, :reference] begin
    using AdditionalDistributions
    using Distributions
    using SpecialFunctions
    using Statistics
    using Test

    @testset "Rademacher" begin
        d = Rademacher()
        @test params(d) == ()
        @test pdf(d, -1) == 0.5
        @test pdf(d, 1) == 0.5
        @test pdf(d, 0) == 0.0
        @test cdf(d, -2) == 0.0
        @test cdf(d, 0) == 0.5
        @test cdf(d, 2) == 1.0
        @test quantile(d, 0.5) == -1
        @test quantile(d, 0.5000001) == 1
        @test mean(d) == 0
        @test var(d) == 1
    end

    @testset "ZIP" begin
        λ, p0 = 2.0, 0.3
        d = ZIP(λ, p0)
        pdf0 = p0 + (1 - p0) * exp(-λ)
        pdf1 = (1 - p0) * exp(-λ) * λ
        pdf2 = (1 - p0) * exp(-λ) * λ^2 / 2

        @test params(d) == (λ, p0)
        @test pdf(d, 0) ≈ pdf0 atol=1e-14
        @test pdf(d, 1) ≈ pdf1 atol=1e-14
        @test cdf(d, 2) ≈ pdf0 + pdf1 + pdf2 atol=1e-14
        @test logpdf(d, 1) ≈ log(pdf1) atol=1e-14
        @test mean(d) ≈ λ * (1 - p0) atol=1e-14
    end

    @testset "ZINB" begin
        r, θ, p = 4, 0.2, 0.7
        d = ZINB(r, θ, p)
        nb = NegativeBinomial(r, p)

        @test params(d) == (r, θ, p)
        @test pdf(d, 0) ≈ θ + (1 - θ) * pdf(nb, 0) atol=1e-14
        @test pdf(d, 3) ≈ (1 - θ) * pdf(nb, 3) atol=1e-14
        @test cdf(d, 3) ≈ θ + (1 - θ) * cdf(nb, 3) atol=1e-14
        @test logpdf(d, 3) ≈ log((1 - θ) * pdf(nb, 3)) atol=1e-14
    end

    @testset "ZIB" begin
        n, θ, p0 = 10, 0.3, 0.25
        d = ZIB(n, θ, p0)
        bin = Binomial(n, θ)

        @test params(d) == (n, θ, p0)
        @test pdf(d, 0) ≈ p0 + (1 - p0) * pdf(bin, 0) atol=1e-14
        @test pdf(d, 3) ≈ (1 - p0) * pdf(bin, 3) atol=1e-14
        @test cdf(d, 3) ≈ p0 + (1 - p0) * cdf(bin, 3) atol=1e-14
    end

    @testset "BetaNegBinomial" begin
        r, α, β = 3, 2.5, 1.25
        d = BetaNegBinomial(r, α, β)
        k = 4
        logp = loggamma(r + k) + logbeta(α + r, β + k) -
               loggamma(k + 1) - loggamma(r) - logbeta(α, β)

        @test params(d) == (r, α, β)
        @test pdf(d, k) ≈ exp(logp) atol=1e-14
        @test logpdf(d, k) ≈ logp atol=1e-14
        @test cdf(d, -1) == 0.0
    end

    @testset "Borel" begin
        a = 0.4
        d = Borel(a)
        pmf(k) = exp(-a * k + (k - 1) * log(a * k) - loggamma(k + 1))

        @test params(d) == (a,)
        @test pdf(d, 3) ≈ pmf(3) atol=1e-14
        @test logpdf(d, 3) ≈ log(pmf(3)) atol=1e-14
        @test cdf(d, 3) ≈ sum(pmf(k) for k in 1:3) atol=1e-14
        @test cdf(d, 0) == 0.0
    end
end
