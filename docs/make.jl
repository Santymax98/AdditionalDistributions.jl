# docs/make.jl
using Pkg
Pkg.develop(path = joinpath(@__DIR__, ".."))
using AdditionalDistributions
using Documenter
using Distributions
using DocumenterVitepress

DocMeta.setdocmeta!(AdditionalDistributions, :DocTestSetup, :(using AdditionalDistributions); recursive=true)

makedocs(;
    modules=[AdditionalDistributions],
    repo = Remotes.GitHub("Santymax98", "AdditionalDistributions.jl"),
    authors="Santiago Jimenez Ramos <santiago.jimenez@ufpe.br>, and contributors",
    sitename="AdditionalDistributions.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo       = "https://github.com/Santymax98/AdditionalDistributions.jl",
        devbranch  = "main",                               # <-- IMPORTANTE
        devurl     = "dev",                                # opcional
        deploy_url = "https://santymax98.github.io/AdditionalDistributions.jl", # opcional
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
    repo         = "github.com/Santymax98/AdditionalDistributions.jl",
    target       = "build",
    devbranch    = "main",
    branch       = "gh-pages",
    push_preview = true,
)