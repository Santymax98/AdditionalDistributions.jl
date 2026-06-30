# docs/make.jl
using Pkg
# desarrolla el paquete localmente (no requiere registro)
Pkg.develop(path = joinpath(@__DIR__, ".."))
using AdditionalDistributions
using Documenter
using GR
using Distributions
using DocumenterVitepress

DocMeta.setdocmeta!(AdditionalDistributions, :DocTestSetup, :(using AdditionalDistributions); recursive=true)

makedocs(;
    modules=[AdditionalDistributions],
    repo = Remotes.GitHub("Santymax98", "AdditionalDistributions.jl"),
    authors="Santiago Jimenez Ramos <santiago.jimenez@ufpe.br>, and contributors",
    sitename="AdditionalDistributions.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/Santymax98/AdditionalDistributions.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Compatibility" => "Compatibility.md",
        "Bestiary" => [
            "Discrete Distributions"=>"bestiary/discrete.md",
            "Continuous Distribution"=>"bestiary/continuous.md",
            "Multivariate Distributions"=>"bestiary/multivariate.md",
        ],
    ],
)


# Minimal sitemap to make the stable documentation easier to discover.
docs_url = "https://santymax98.github.io/AdditionalDistributions.jl/stable"
docs_pages = [
    "$(docs_url)/",
    "$(docs_url)/getting_started",
    "$(docs_url)/Compatibility",
    "$(docs_url)/bestiary/discrete",
    "$(docs_url)/bestiary/continuous",
    "$(docs_url)/bestiary/multivariate",
]

open(joinpath(@__DIR__, "build", "sitemap.xml"), "w") do io
    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    println(io, "<urlset xmlns=\"http://www.sitemaps.org/schemas/sitemap/0.9\">")
    for url in docs_pages
        println(io, "  <url>")
        println(io, "    <loc>$(url)</loc>")
        println(io, "  </url>")
    end
    println(io, "</urlset>")
end

DocumenterVitepress.deploydocs(;
    repo = "github.com/Santymax98/AdditionalDistributions.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch="main",
    branch = "gh-pages",
    push_preview = true,
)
