@testitem "Student radial exact fast paths" begin
    using Test
    using Distributions
    using AdditionalDistributions

    const AD = AdditionalDistributions

    # ---------------------------------------------------------
    # Ordinary range:
    # compare against the generic Distributions quantile.
    # ---------------------------------------------------------
    ps_regular = Float64[
        1e-12,
        1e-10,
        1e-8,
        1e-6,
        1e-4,
        0.01,
        0.1,
        0.5,
        0.9,
        0.99,
        1 - 1e-8,
        1 - 1e-12,
    ]

    for ν in (1.0, 2.0, 4.0)
        χ = Chisq(ν)
        invsqrtν = inv(sqrt(ν))

        for p in ps_regular
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

    # ---------------------------------------------------------
    # ν = 1:
    #
    # χ²₁ = Z², hence for the radial scale
    #
    #   R = sqrt(Qχ²₁(p))
    #     = sqrt(2) * erfinv(p).
    #
    # Do NOT use Chisq.quantile as the reference in the
    # extreme lower tail: its relative accuracy there is not
    # suitable for validating this analytical fast path.
    #
    # Instead verify the inverse identity
    #
    #   p = erf(R / sqrt(2)).
    # ---------------------------------------------------------
    χ1 = Chisq(1.0)

    ps_nu1_tail = Float64[
        1e-300,
        1e-200,
        1e-100,
        1e-50,
        1e-30,
        1e-20,
        1e-16,
        1e-14,
    ]

    let previous = 0.0
        for p in ps_nu1_tail
            r = AD._scale_t(χ1, 1.0, p)

            @test isfinite(r)
            @test r >= 0
            @test r >= previous

            p_back = AD.SpecialFunctions.erf(r / sqrt(2.0))

            @test isapprox(
                p_back,
                p;
                rtol=5e-13,
                atol=0.0,
            )

            previous = r
        end
    end

    # ---------------------------------------------------------
    # ν = 2:
    #
    # χ²₂ ~ Exponential(scale=2), therefore
    #
    #   sqrt(Qχ²₂(p) / 2) = sqrt(-log(1-p)).
    # ---------------------------------------------------------
    χ2 = Chisq(2.0)

    for p in ps_regular
        fast = AD._scale_t(
            χ2,
            inv(sqrt(2.0)),
            p,
        )

        exact = sqrt(-log1p(-p))

        @test fast == exact
    end

    # ---------------------------------------------------------
    # ν = 4:
    #
    # Verify the Lambert-W implementation throughout the
    # ordinary range and deep lower/upper tails.
    #
    # The smallest Float64 value deliberately exercises the
    # generic χ² fallback.
    # ---------------------------------------------------------
    χ4 = Chisq(4.0)

    ps_nu4 = Float64[
        nextfloat(0.0),
        1e-300,
        1e-200,
        1e-100,
        1e-50,
        1e-30,
        1e-20,
        1e-16,
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

    for p in ps_nu4
        fast = AD._scale_t(χ4, 0.5, p)
        ref = sqrt(quantile(χ4, p) / 4.0)

        @test isfinite(fast)
        @test fast >= 0

        @test isapprox(
            fast,
            ref;
            rtol=5e-13,
            atol=1e-14,
        )
    end
end

@testitem "Student-t radial scaling keeps generic Real dispatch" begin
    using AdditionalDistributions
    using Distributions
    using ForwardDiff

    const AD = AdditionalDistributions

    u = ForwardDiff.Dual{Nothing}(0.37, 1.0)
    invsqrtν = one(u) / sqrt(5.0)

    @test applicable(AD._scale_t, 5.0, u)
    @test applicable(AD._scale_t, Chisq(5.0), invsqrtν, u)
end
