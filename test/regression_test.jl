@testitem "Regression – corrected discrete semantics" tags=[:regression, :discrete] begin
    using Distributions: params, pdf, logpdf, cdf

    d = Rademacher()
    @test params(d) == ()
    @test pdf(d, 0) == 0.0
    @test logpdf(d, 0) == -Inf
    @test cdf(d, -2) == 0.0
    @test cdf(d, 0) == 0.5
    @test cdf(d, 2) == 1.0

    z = Zeta(2.0)
    @test params(z) == (2.0,)
    @test pdf(z, 1.5) == 0.0
    @test logpdf(z, 1.5) == -Inf
    @test cdf(z, Inf) == 1.0

    zinb = ZINB(4, 0.2, 0.7)
    @test zinb isa ZINB
    @test params(zinb) == (4, 0.2, 0.7)
    @test pdf(zinb, 0) ≈ 0.2 + 0.8 * 0.7^4
    @test pdf(zinb, 1.2) == 0.0
end

@testitem "Regression – Delaporte normalization" tags=[:regression, :discrete] begin
    using Distributions: pdf, logpdf

    d = Delaporte(2.5, 2.0, 1.0)
    s = sum(pdf(d, k) for k in 0:120)

    @test isapprox(s, 1.0; atol=1e-8)
    @test isapprox(logpdf(d, 4), log(pdf(d, 4)); atol=1e-12)
    @test logpdf(d, -1) == -Inf
end

@testitem "Regression – continuous corrections" tags=[:regression, :continuous] begin
    using Distributions: cdf, quantile
    using Statistics: var

    bs = BirnbaumSaunders(0.0, 1.0, 1.0)
    @test isfinite(var(bs))

    p = PERT(0.0, 1.0, 3.0)
    @test cdf(p, -1.0) == 0.0
    @test cdf(p, 4.0) == 1.0
    @test quantile(p, 0.0) == 0.0
    @test quantile(p, 1.0) == 3.0

    cb = CrystalBall(1.5, 3.5, 0.0, 1.0)
    @test cdf(cb, -Inf) == 0.0
    @test cdf(cb, Inf) == 1.0
    @test cdf(cb, quantile(cb, 0.25)) ≈ 0.25 atol=1e-8
    @test cdf(cb, quantile(cb, 0.75)) ≈ 0.75 atol=1e-8
end