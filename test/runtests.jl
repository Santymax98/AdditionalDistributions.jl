using TestItems
using TestItemRunner
# --- OPCIÓN A: correr TODO ---
@run_package_tests

# --- OPCIÓN B: solo las pruebas con tag :Borel ---
#@run_package_tests (filter = ti -> :Alpha in ti.tags)
