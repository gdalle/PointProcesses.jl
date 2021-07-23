"""
    PointProcess{L, M}

The common supertype of all point processes with location type `L` and mark type `M`.
"""
abstract type PointProcess{L,M} end

"""
    params(pp)

Retrieve point process parameters as a `NamedTuple`.
"""
StatsBase.params(pp::PointProcess) = ntfromstruct(pp)
