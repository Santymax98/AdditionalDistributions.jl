"""
    CDFResult(value, error, inform, neval, algorithm)

Structured result returned by `cdf_result` for multivariate rectangular
probabilities.

Fields:
- `value`: estimated probability.
- `error`: estimated absolute integration error.
- `inform`: convergence/status code.
- `neval`: approximate number of integrand evaluations.
- `algorithm`: integration algorithm identifier.
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