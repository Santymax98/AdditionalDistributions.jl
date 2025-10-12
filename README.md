# AdditionalDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl)


The [*AdditionalDistributions*](https://github.com/Santymax98/AdditionalDistributions.jl) package is a comprehensive extension of [*Distributions*](https://github.com/JuliaStats/Distributions.jl). It is designed to expand the functionality of the base package by incorporating both discrete and continuous probability distributions that are not included in `Distributions`, either due to their specialized nature or because they are less commonly used. *AdditionalDistributions* aims to provide a broader range of statistical tools for data analysis, simulation, and probabilistic modeling, catering to both academic and scientific needs where these additional distributions are essential.

**Purpose and Scope**

The purpose of `AdditionalDistributions.jl` is to serve as a well-maintained repository for distributions that are too exotic or specialized to be included in the main `Distributions.jl` package. This includes distributions that are frequently required in specific fields of research but are not yet available in Julia, as well as new or experimental distributions that may not be widely used but have significant potential applications.

For instance, you can define the following distributions (among many others):

```julia
julia> using Distributions, AdditionalDistributions
julia> BetaNegBinomial(r, α, β) # Discrete univariate
julia> Lomax(α, λ)              # Continuous univariate
julia> ZINB(r, θ, p)            # Discrete univariate
julia> Gompertz(η, b)           # Continuous univariate
```

These distributions, along with others like the ARGUS or Zero-Inflated Poisson (ZIP), are essential in various fields of research but are not included in the base `Distributions.jl` package. `AdditionalDistributions.jl` provides these and more, ensuring that users have access to a broad spectrum of statistical tools.

**Key Features**

- **Extensive Range of Distributions**: The package includes a wide variety of distributions, some of which are not commonly found in standard statistical libraries.
- **Seamless Integration**: All distributions in `AdditionalDistributions.jl` are fully compatible with `Distributions.jl`, making it easy to use them in conjunction with other packages in the Julia ecosystem.
- **Detailed Documentation**: We are committed to providing thorough documentation, including usage examples, theoretical background, and references for each distribution. This will ensure that users can effectively implement and understand the distributions available in the package.
- **Community-Driven Development**: We welcome and encourage community contributions, whether it's proposing new distributions, improving documentation, or contributing code to further expand the package.

**Maintaining Relevance and Utility**

Given the importance of these additional distributions in various academic and scientific applications, we are committed to the ongoing maintenance and development of `AdditionalDistributions.jl`. Our goal is to ensure that the package remains a valuable resource for the Julia community. To support this, we are considering moving the package under the `JuliaStats` organization, which would provide additional support and visibility, ensuring its long-term sustainability.

With *AdditionalDistributions*, you can:

- **Sample from distributions:** Draw random samples from a variety of distributions.
- **Calculate moments and other properties:** Obtain moments (such as mean, variance, skewness, and kurtosis), entropy, and other statistical properties.
- **Evaluate probability density/mass functions:** Compute the probability density functions (pdf) and their logarithms (logpdf).
- **Utilize moment-generating, quantile, and characteristic functions:** Access moment-generating functions, quantile functions, and characteristic functions for in-depth statistical analysis.

**Future Directions**

In the future, we plan to implement maximum likelihood estimators and potentially introduce additional multivariate distributions to further enrich the package.