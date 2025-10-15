# Compatibility

## Compatibility with Distributions.jl

The *ProbabilityDistributions* package is fully compatible with the [*Distributions.jl*](https://github.com/JuliaStats/Distributions.jl) package. As an extension of *Distributions.jl*, *ProbabilityDistributions* supports all the core functionalities provided by *Distributions.jl*, including but not limited to:

### Truncated Distributions 

Generate and work with truncated versions of existing distributions.

```julia
julia> d = Alpha()
Alpha{Float64}(α=1.0, β=1.0)
julia> d_truncated = Truncated(d, 0.0, 1.0)
Truncated(Alpha{Float64}(α=1.0, β=1.0); lower=0.0, upper=1.0)
```
### Censored Distributions

Handle censored data and perform analyses accordingly.

```julia
julia> d1 = Gompertz()
Gompertz{Float64}(η=1.0, b=1.0)
d_censored = censored(d1, 1.5, 10)
Censored(Gompertz{Float64}(η=1.0, b=1.0); lower=1.5, upper=10.0)
```

### Mixture Distributions

Create and analyze mixtures of different probability distributions.

```julia
julia> d_mixture = MixtureModel(Maxwell, [2.0, 1.0, 5.5], [0.2, 0.5, 0.3])
MixtureModel{Maxwell}(K = 3)
components[1] (prior = 0.2000): Maxwell{Float64}(a=2.0)
components[2] (prior = 0.5000): Maxwell{Float64}(a=1.0)
components[3] (prior = 0.3000): Maxwell{Float64}(a=5.5)
```

### Order Statistics

Compute and work with order statistics.
```julia
julia> OrderStatistic(Burr(), 10, 1)
OrderStatistic{Burr{Float64}, Continuous}(
dist: Burr{Float64}(c=1.0, k=1.0, λ=1.0)
n: 10
rank: 1
)

julia> OrderStatistic(Logarithmic(), 10, 5)
OrderStatistic{Logarithmic{Float64}, Discrete}
(dist: Logarithmic{Float64}(a=0.5) 
n: 10 
rank: 5
)
```


## Integration with Other Packages

Because of its compatibility with *Distributions.jl*, *ProbabilityDistributions* seamlessly integrates with other Julia packages that also build on *Distributions.jl*. This includes:

- **[Turing.jl](https://github.com/TuringLang/Turing.jl)**: For probabilistic programming and Bayesian inference.
- **[HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl)**: For hypothesis testing and statistical testing.
- **[Copulas.jl](https://github.com/lrnv/Copulas.jl)**: For copula-based modeling and simulations.

This broad compatibility ensures that you can use *ProbabilityDistributions* in a wide range of statistical applications, from advanced simulations to Bayesian analysis and hypothesis testing.