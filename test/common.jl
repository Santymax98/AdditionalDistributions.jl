@testmodule M begin
    using Test
    using Random, StableRNGs
    using Distributions
    using StatsBase
    using AdditionalDistributions

    const rng = StableRNG(2024)

    is_discrete(d)   = d isa Distributions.DiscreteUnivariateDistribution
    is_continuous(d) = d isa Distributions.ContinuousUnivariateDistribution

    _minmax(d) = (try Distributions.minimum(d) catch; -Inf end, try Distributions.maximum(d) catch;  Inf end)

    function check_cdf_edges(d::Distributions.UnivariateDistribution; atol=1e-12)
        @info "Testing CDF edges for $d..."

        # Universal CDF/CCDF limits. These should hold independently of the support.
        @test Distributions.cdf(d, -Inf) ≈ 0.0 atol=atol
        @test Distributions.cdf(d,  Inf) ≈ 1.0 atol=atol
        @test Distributions.ccdf(d, -Inf) ≈ 1.0 atol=atol
        @test Distributions.ccdf(d,  Inf) ≈ 0.0 atol=atol

        # NaN should propagate instead of being silently converted into 0 or 1.
        @test isnan(Distributions.cdf(d, NaN))
        @test isnan(Distributions.ccdf(d, NaN))

        # Use generic probes plus support-dependent probes.
        a, b = _minmax(d)
        xs = Float64[-Inf, -10.0, -1.0, -0.0, 0.0, 1e-12, 1.0, 10.0, Inf]

        if isfinite(a)
            push!(xs, prevfloat(float(a)), float(a), nextfloat(float(a)))
        end
        if isfinite(b)
            push!(xs, prevfloat(float(b)), float(b), nextfloat(float(b)))
        end

        sort!(unique!(xs))
        vals = Distributions.cdf.(Ref(d), xs)

        for v in vals
            if !isnan(v)
                @test -atol <= v <= 1.0 + atol
            end
        end

        finite_vals = filter(!isnan, vals)
        @test all(diff(finite_vals) .>= -atol)

        # -0.0 and 0.0 represent the same real number.
        @test Distributions.cdf(d, -0.0) ≈ Distributions.cdf(d, 0.0) atol=atol

        if isfinite(a) && a == 0
            @test Distributions.cdf(d, -eps(Float64)) ≈ 0.0 atol=atol
            @test -atol <= Distributions.cdf(d, -0.0) <= 1.0 + atol
        end

        nothing
    end

    function check_discrete(d; nsample=2_000, nprobe=50)
        @info "Testing $d..."
        @test d isa Distributions.DiscreteUnivariateDistribution

        check_cdf_edges(d)

        X = rand(rng, d, nsample)

        @test all(x -> Distributions.insupport(d, x), X)

        cdfX = Distributions.cdf.(Ref(d), X)
        @test all(y -> 0.0 <= y <= 1.0, cdfX)

        xs  = sort!(unique(X))
        xs  = xs[1:min(end, nprobe)]
        cfs = Distributions.cdf.(Ref(d), xs)
        @test issorted(cfs)

        pick = StatsBase.sample(X, min(nprobe, length(X)); replace=false)

        logp = Distributions.logpdf.(Ref(d), pick)
        p    = Distributions.pdf.(Ref(d), pick)

        @test all(isfinite, logp)
        @test all(p .>= 0)
        @test all(abs.(exp.(logp) .- p) .<= 1e-12)

        F  = Distributions.cdf.(Ref(d), pick)
        pm = Distributions.pdf.(Ref(d), pick)
        ps = (F .+ (F .- pm)) ./ 2

        @test Distributions.quantile.(Ref(d), ps) == pick

        a, b = _minmax(d)

        if isfinite(a)
            a_below = a isa Integer ? a - one(a) : prevfloat(float(a))
            @test Distributions.cdf(d, a_below) ≈ 0.0 atol=1e-12
            @test Distributions.cdf(d, a) ≈ Distributions.pdf(d, a) atol=1e-12
        end

        if isfinite(b)
            @test Distributions.cdf(d, b) ≈ 1.0 atol=1e-12
        else
            q = Distributions.quantile(d, 0.999)
            @test Distributions.cdf(d, q) ≥ 0.999 - 1e-12
        end

        nothing
    end

    function check_continuous(d; nsample=2_000, nprobe=50)
        @info "Testing $d..."
        @test d isa Distributions.ContinuousUnivariateDistribution

        check_cdf_edges(d)

        X = rand(rng, d, nsample)

        @test all(isfinite, X)
        @test all(x -> Distributions.insupport(d, x), X)

        cdfX = Distributions.cdf.(Ref(d), X)
        @test all(y -> 0.0 <= y <= 1.0, cdfX)

        xs  = sort!(copy(X))
        xs  = xs[1:min(end, nprobe)]
        cfs = Distributions.cdf.(Ref(d), xs)
        @test issorted(cfs)

        pick = StatsBase.sample(X, min(nprobe, length(X)); replace=false)

        logp = Distributions.logpdf.(Ref(d), pick)
        p    = Distributions.pdf.(Ref(d), pick)

        @test all(isfinite, logp)
        @test all(p .>= 0)
        @test all(abs.(exp.(logp) .- p) .<= 1e-10)

        ps = Distributions.cdf.(Ref(d), pick)
        ps = clamp.(ps, 1e-12, 1 - 1e-12)
        qs = Distributions.quantile.(Ref(d), ps)

        @test all(isfinite, qs)

        ps_roundtrip = Distributions.cdf.(Ref(d), qs)
        tol = sqrt(eps(Float64))
        @test all(isapprox(p2, p; atol=tol, rtol=tol) for (p2, p) in zip(ps_roundtrip, ps))

        a, b = _minmax(d)
        if isfinite(a)
            @test Distributions.cdf(d, a) ≈ 0.0 atol=1e-12
        end
        if isfinite(b)
            @test Distributions.cdf(d, b) ≈ 1.0 atol=1e-12
        end

        nothing
    end

    function check(d::Distributions.UnivariateDistribution; kwargs...)
        if is_discrete(d)
            check_discrete(d; kwargs...)
        else
            check_continuous(d; kwargs...)
        end
    end
end