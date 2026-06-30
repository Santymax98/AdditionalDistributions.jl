# Getting Started

## Installation

The *AdditionalDistributions* package is available through the Julia package system by running `Pkg.add("AdditionalDistributions")`.
Throughout, we assume that you have installed the package.

## Starting With a Poisson Zero-inflated Distribution ZIP

We start by drawing 100 observations from a Poisson zero-inflated with parameters `λ = 5.0` and `p = 0.2` random variable.

The first step is to set up the environment:

```julia
julia> using Random, Distributions, AdditionalDistributions

julia> Random.seed!(2024) # Setting the seed
```

Then, we create a Poisson zero-inflated distribution `d` and obtain samples using `rand`:

```julia
julia> d = ZIP(5.0, 0.2)
ZIP{Float64}(λ=5.0, p=0.2)
```

The object `d` represents a probability distribution, in our case the Poisson zero-inflated distribution.
One can query its properties such as the mean:

```julia
julia> mean(d)
4.0
```

We can also draw samples from `d` with `rand`.
```julia
julia> samples = rand(d, 100)
100-element Vector{Int64}:
 0
 3
 6
 4
 5
 ⋮
```

You can easily obtain the `pdf`, `cdf`, `quantile`, and many other functions for a distribution. For instance, the median (50th percentile) and the 95th percentile for the Poisson zero inflated distribution are given by:

```julia
julia> quantile.(ZIP(), [0.5, 0.95])
2-element Vector{Int64}:
 0
 2
```

## Using Other Distributions

The package contains a large number of `discrete` and `continuous` distributions in addition to those implemented in Distributions.jl.

For instance, you can define the following distributions (among many others):

```julia
julia> BetaNegBinomial(r, α, β) # Discrete univariate
julia> Lomax(α, λ)              # Continuous univariate
julia> ZINB(r, θ, p)            # Discrete univariate
julia> Gompertz(η, b)           # Continuous univariate
```