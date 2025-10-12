@testset "Aqua.jl" begin
    Aqua.test_all(
        AdditionalDistributions;
        persistent_tasks = VERSION != v"1.10.10", # workaround for Julia 1.10.10
        ambiguities = false,
    )
end

