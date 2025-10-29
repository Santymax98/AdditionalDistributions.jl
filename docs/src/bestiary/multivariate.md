## Extra Multivariate Distributions

The following multivariate distributions extend the behavior of `Distributions.jl` to include **numerically stable cumulative
probability computations** using custom Quasi–Monte Carlo (QMC) integration routines.

---

```@docs
MvGaussian
```

**Notes**

* Equivalent to `Distributions.MvNormal`, but implements its own `cdf(a, b)` method based on a high–precision QMC integrator.
* Compatible with reference datasets from [MvNormalCDF.jl](https://github.com/JuliaStats/MvNormalCDF.jl), adjusted for local reproducibility.
* Fully compatible with `mean`, `cov`, `pdf`, `logpdf`, and `rand`.

**Example**

```julia
μ = [0.0, 0.0]
Σ = [1.0 0.5; 0.5 1.0]
d = MvGaussian(μ, Σ)
cdf(d, [-1.0, -1.0], [1.0, 1.0])
```

```@docs
MvTStudent
```

**Notes**

* Equivalent to `Distributions.MvTDist`, but with a custom `cdf(a, b)` based on the same QMC integrator.
* For ν → ∞, numerical results converge to those of `MvGaussian`.
* Benchmarked against `pmvt` from *mvtnorm* (R, Genz & Bretz, 2002).

**Example**

```julia
ν = 10
Σ = [1.0 0.4; 0.4 1.0]
d = MvTStudent(ν, Σ)
cdf(d, [-1.0, -1.0], [1.0, 1.0])
```

---

### Implementation Notes

Both types serve as lightweight wrappers around their `Distributions.jl` counterparts:

```julia
struct MvGaussian{D<:Distributions.MvNormal}
    dist::D
end

struct MvTStudent{D<:Distributions.MvTDist}
    dist::D
end
```

All statistical methods (`mean`, `cov`, `pdf`, etc.) are delegated to their internal distribution objects, ensuring full compatibility while
adding numerical integration support for the cumulative distribution function.

```@docs
AdditionalDistributions.mvtcdf
```