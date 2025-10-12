@testmodule M begin
    using Test
    using Random, StableRNGs
    using Distributions
    using StatsBase
    using AdditionalDistributions

    const rng = StableRNG(2024)

    is_discrete(d)   = d isa Distributions.DiscreteUnivariateDistribution
    is_continuous(d) = d isa Distributions.ContinuousUnivariateDistribution

    _minmax(d) = (try Distributions.minimum(d) catch; -Inf end,
                  try Distributions.maximum(d) catch;  Inf end)

    function check_discrete(d; nsample=2_000, nprobe=50)
        @info "Testing $d..."
        @test d isa Distributions.DiscreteUnivariateDistribution

        X = rand(rng, d, nsample)

        # Usa insupport() (no "x ∈ support(d)") para evitar problemas con Inf
        @test all(x -> Distributions.insupport(d, x), X)

        # CDF en [0,1] y monótona
        @test all(0.0 .<= Distributions.cdf.(Ref(d), X) .<= 1.0)
        xs  = sort!(unique(X))
        xs  = xs[1:min(end, nprobe)]
        cfs = Distributions.cdf.(Ref(d), xs)
        @test issorted(cfs)

        # logpdf/pdf consistentes y pmf ≥ 0
        pick = StatsBase.sample(X, min(nprobe, length(X)); replace=false)
        @test all(isfinite.(Distributions.logpdf.(Ref(d), pick)))
        @test all(Distributions.pdf.(Ref(d), pick) .>= 0)
        @test all(abs.(exp.(Distributions.logpdf.(Ref(d), pick)) .- Distributions.pdf.(Ref(d), pick)) .<= 1e-12)

        # Round-trip de cuantiles usando punto medio del salto
        F  = Distributions.cdf.(Ref(d), pick)
        pm = Distributions.pdf.(Ref(d), pick)
        ps = (F .+ (F .- pm)) ./ 2
        @test Distributions.quantile.(Ref(d), ps) == pick

        # Bordes (mín/max) si son finitos
        a, b = _minmax(d)
        if isfinite(a)
            a_below = isinteger(a) ? Int(a) - 1 : prevfloat(a)
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

        X = rand(rng, d, nsample)

        @test all(isfinite, X)
        @test all(x -> Distributions.insupport(d, x), X)
        @test all(0.0 .<= Distributions.cdf.(Ref(d), X) .<= 1.0)

        xs  = sort!(copy(X))
        xs  = xs[1:min(end, nprobe)]
        cfs = Distributions.cdf.(Ref(d), xs)
        @test issorted(cfs)

        pick = StatsBase.sample(X, min(nprobe, length(X)); replace=false)
        @test all(isfinite.(Distributions.logpdf.(Ref(d), pick)))
        @test all(abs.(exp.(Distributions.logpdf.(Ref(d), pick)) .- Distributions.pdf.(Ref(d), pick)) .<= 1e-10)

        ps = Distributions.cdf.(Ref(d), pick)
        ps = clamp.(ps, 1e-12, 1 - 1e-12)
        qs = Distributions.quantile.(Ref(d), ps)
        @test all(isfinite, qs)
        @test maximum(abs.(qs .- pick)) ≤ 1e-6

        a, b = _minmax(d)
        if isfinite(a); @test Distributions.cdf(d, a) ≈ 0.0 atol=1e-12; end
        if isfinite(b); @test Distributions.cdf(d, b) ≈ 1.0 atol=1e-12; end
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