```@setup plotdensity
using Distributions, GR

# display figures as SVGs
GR.inline("svg")

# plot probability density of continuous distributions
function plotdensity(
    (xmin, xmax),
    dist::ContinuousUnivariateDistribution;
    npoints=299,
    title="",
    kwargs...,
)
    figure(;
        title=title,
        xlabel="x",
        ylabel="density",
        grid=false,
        backgroundcolor=0, # white instead of transparent background for dark Documenter scheme
        font="Helvetica_Regular", # work around https://github.com/JuliaPlots/Plots.jl/issues/2596
        linewidth=2.0, # thick lines
        kwargs...,
    )
    return plot(range(xmin, xmax; length=npoints), Base.Fix1(pdf, dist))
end

# convenience function with automatic title
function plotdensity(
    xmin_xmax,
    ::Type{T},
    args=();
    title=string(T) * "(" * join(args, ", ") * ")",
    kwargs...
) where {T<:ContinuousUnivariateDistribution}
    return plotdensity(xmin_xmax, T(args...); title=title, kwargs...)
end
```
## Extra Continuous Distributions
```@docs
Alpha
```
```@example plotdensity
using Distributions, AdditionalDistributions
plotdensity((0.0, 1.0), Alpha, (1, 1)) # hide
```
```@docs
Argus
```
```@example plotdensity
plotdensity((0.0, 1.0), Argus, (1, 1)) # hide
```
```@docs
Benini
```
```@example plotdensity
plotdensity((1.0, 2.0), Benini, (1, 1, 1)) # hide
```
```@docs
Benktander_Type1
```
```@example plotdensity
plotdensity((1.0, 2.0), Benktander_Type1, (1, 1)) # hide
```
```@docs
Benktander_Type2
```
```@example plotdensity
plotdensity((1.0, 2.0), Benktander_Type2, (1, 1)) # hide
```
```@docs
Bhattacharjee
```
```@example plotdensity
plotdensity((0.0, 1.0), Bhattacharjee, (0, 1, 1)) # hide
```
```@docs
BirnbaumSaunders
```
```@example plotdensity
plotdensity((0.0, 1.0), BirnbaumSaunders, (0, 1, 1)) # hide
```
```@docs
Bradford
```
```@example plotdensity
plotdensity((0.0, 1.0), Bradford, (0.5,)) # hide
```
```@docs
Burr
```
```@example plotdensity
plotdensity((0.0, 1.0), Burr, (1,1,1)) # hide
```
```@docs
CrystalBall
```
```@example plotdensity
plotdensity((0.0, 1.0), CrystalBall, (1, 2 ,0, 1)) # hide
```
```@docs
Dagum
```
```@example plotdensity
plotdensity((0.0, 1.0), Dagum, (1,1,1)) # hide
```
```@docs
Gompertz
```
```@example plotdensity
plotdensity((0.0, 1.0), Gompertz, (1,1)) # hide
```
```@docs
IrwinHall
```
```@example plotdensity
plotdensity((0.0, 3.0), IrwinHall, (3)) # hide
```
```@docs
Lomax
```
```@example plotdensity
plotdensity((0.0, 1.0), Lomax, (1,1)) # hide
```
```@docs
Maxwell
```
```@example plotdensity
plotdensity((0.0, 1.0), Maxwell, (1,)) # hide
```
```@docs
Nakagami
```
```@example plotdensity
plotdensity((0.0, 1.0), Nakagami, (0.5, 1.0)) # hide
```
```@docs
PERT
```
```@example plotdensity
plotdensity((0.0, 1.0), PERT, (0.0, 0.5, 1.0)) # hide
```
```@index
Pages = ["Distributions.md"]
```