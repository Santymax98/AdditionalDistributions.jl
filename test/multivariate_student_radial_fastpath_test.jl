@testitem "Student radial exact fast paths" begin
    using Test
    using Distributions
    using AdditionalDistributions

    const AD = AdditionalDistributions

    ps = [
        eps(Float64),
        1e-12,
        1e-8,
        1e-4,
        0.01,
        0.1,
        0.5,
        0.9,
        0.99,
        1 - 1e-8,
        1 - 1e-12,
        prevfloat(1.0),
    ]

    for ν in (1.0, 2.0)
        χ = Chisq(ν)
        invsqrtν = inv(sqrt(ν))

        for p in ps
            fast = AD._scale_t(χ, invsqrtν, p)
            ref = sqrt(quantile(χ, p) / ν)

            @test isfinite(fast)
            @test fast >= 0
            @test isapprox(
                fast,
                ref;
                rtol=1e-11,
                atol=1e-14,
            )
        end
    end
end
