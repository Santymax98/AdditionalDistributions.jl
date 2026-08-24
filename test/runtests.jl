using Test
using TestItems
using TestItemRunner
using Aqua
using AdditionalDistributions

@testset "Aqua.jl" begin
    Aqua.test_all(
        AdditionalDistributions;
        persistent_tasks = VERSION != v"1.10.10", # workaround Julia 1.10.10 bug
        ambiguities = false,
    )
end

# Run the complete TestItems suite by default.
@run_package_tests #(filter = ti -> :discrete in ti.tags)

# Useful local filters for development:
# @run_package_tests (filter = ti -> :multivariate in ti.tags)
# @run_package_tests (filter = ti -> :continuous in ti.tags)
# @run_package_tests (filter = ti -> :discrete in ti.tags)
# @run_package_tests (filter = ti -> :reference in ti.tags)
# @run_package_tests (filter = ti -> !(:slow in ti.tags))
