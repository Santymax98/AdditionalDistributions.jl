using Test
using TestItems
using TestItemRunner
using Aqua
using AdditionalDistributions

# --- Pruebas Aqua.jl ---
@testset "Aqua.jl" begin
    Aqua.test_all(
        AdditionalDistributions;
        persistent_tasks = VERSION != v"1.10.10", # workaround Julia 1.10.10 bug
        ambiguities = false,  # opcional: puedes poner true si estás depurando métodos
    )
end

# --- Ejecución de TestItems ---

# OPCIÓN A: correr TODO el paquete
@run_package_tests

# --- OPCIÓN B (descomenta para filtrar por tag específico) ---
# @run_package_tests (filter = ti -> :Alpha in ti.tags)
# @run_package_tests (filter = ti -> :Borel in ti.tags)
