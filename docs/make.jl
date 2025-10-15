# docs/make.jl
using Pkg
# desarrolla el paquete localmente (no requiere registro)
Pkg.develop(path = joinpath(@__DIR__, ".."))
using ProbabilityDistributions
using Documenter
using GR
using Distributions
using DocumenterVitepress

DocMeta.setdocmeta!(ProbabilityDistributions, :DocTestSetup, :(using ProbabilityDistributions); recursive=true)

makedocs(;
    modules=[ProbabilityDistributions],
    repo = Remotes.GitHub("Santymax98", "ProbabilityDistributions.jl"),
    authors="Santiago Jimenez Ramos <santiago.jimenez@ufpe.br>, and contributors",
    sitename="ProbabilityDistributions.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/Santymax98/ProbabilityDistributions.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Compatibility" => "Compatibility.md",
        "Bestiary" => [
            "Discrete Distributions"=>"bestiary/discrete.md",
            "Continuous Distribution"=>"bestiary/continuous.md",
        ],
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/Santymax98/ProbabilityDistributions.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch="main",
        branch = "gh-pages",
    push_preview = true,
)
