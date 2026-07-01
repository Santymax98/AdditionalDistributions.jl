"""
    CDFResult(value, error, inform, neval, algorithm)

Structured result returned by [`cdf_result`](@ref) for multivariate rectangular
probabilities.

# Fields

- `value`: estimated probability.
- `error`: estimated absolute integration error.
- `inform`: convergence/status code.
- `neval`: requested integration budget.
- `algorithm`: integration algorithm identifier.

# `inform` codes

- `0`: estimated error is within tolerance.
- `1`: estimated error is above tolerance for the current budget.
- `2`: invalid dimension.
- `3`: matrix appears not positive semidefinite during preparation.

`CDFResult` can be destructured as `(value, error, inform)` for compatibility
with the legacy `full=true` tuple output.
"""
struct CDFResult{T<:Real}
    value::T
    error::T
    inform::Int
    neval::Int
    algorithm::Symbol
end

Base.iterate(r::CDFResult) = (r.value, 2)

function Base.iterate(r::CDFResult, state::Int)
    state == 2 && return (r.error, 3)
    state == 3 && return (r.inform, 4)
    return nothing
end

Base.length(::CDFResult) = 3

Base.getindex(r::CDFResult, i::Int) =
    i == 1 ? r.value :
    i == 2 ? r.error :
    i == 3 ? r.inform :
    throw(BoundsError(r, i))

Base.first(r::CDFResult) = r.value
