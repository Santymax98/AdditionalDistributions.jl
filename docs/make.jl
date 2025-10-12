push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using AdditionalDistributions
using Documenter
using GR
using Distributions

DocMeta.setdocmeta!(AdditionalDistributions, :DocTestSetup, :(using AdditionalDistributions); recursive=true)

makedocs(;
    modules=[AdditionalDistributions],
    authors="Santiago Jimenez Ramos",
    sitename="AdditionalDistributions.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://Santymax98.github.io/AdditionalDistributions.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "index.md",
        "getting_started.md",
        "Compatibility.md",
        "Distributions.md"
    ],
)

deploydocs(;
    repo="github.com/Santymax98/AdditionalDistributions.jl",
    devbranch="main",
)
