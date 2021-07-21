"""
    PointProcess{M}

The common supertype of all point processes with mark type `M`.
"""
abstract type PointProcess{M} end

"""
    Parameter

An alias for `ComponentVector`, which is how we store point process parameters `θ`.
"""
const Parameter = ComponentVector

Base.eltype(::Type{<:PointProcess{M}}) where {M} = History{M}

StatsBase.params(pp::PointProcess) = ComponentVector(ntfromstruct(pp))
