````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: AdditionalDistributions.jl
  text:
  tagline: An extension of Distributions.jl package. 
  image:
    src: logo.png
    alt: AdditionalDistributions.jl
  actions:
    - theme: brand
      text: Getting started
      link: getting_started
    - theme: alt
      text: View on Github
      link: https://github.com/Santymax98/AdditionalDistributions.jl
    - theme: alt
      text: Bestiary
      link: /bestiary/continuous
---
````
# Welcome to AdditionalDistributions.jl

The [*AdditionalDistributions.jl*](https://github.com/Santymax98/AdditionalDistributions.jl) package is a comprehensive extension of [*Distributions.jl*](https://github.com/JuliaStats/Distributions.jl), designed to expand its functionality by incorporating both discrete and continuous probability distributions not available in the base package. *AdditionalDistributions* aims to enhance the breadth of statistical tools available for data analysis, simulation, and probabilistic modeling.

With *AdditionalDistributions*, you can:

- **Sample from distributions:** Draw random samples from a variety of distributions.
- **Calculate moments and other properties:** Obtain moments (such as mean, variance, skewness, and kurtosis), entropy, and other statistical properties.
- **Evaluate probability density/mass functions:** Compute the probability density functions (pdf) and their logarithms (logpdf).
- **Utilize moment-generating, quantile, and characteristic functions:** Access moment-generating functions, quantile functions, and characteristic functions for in-depth statistical analysis.

In the future, we plan to implement maximum likelihood estimators and potentially introduce additional multivariate distributions to further enrich the package. Collaborators are also welcome.