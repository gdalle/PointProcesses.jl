"""
    PointProcess{M}

The common supertype of all point processes with mark type `M`.
"""
abstract type PointProcess{M} end

"""
    TimedPointProcess{M,P<:PointProcess{M},R<:Real}

Abstract point process with built-in start and end times.

# Fields
- `pp::P`: underlying point process
- `tmin::R`: start time
- `tmax::R`: end time
"""
struct TimedPointProcess{M,P<:PointProcess{M},R<:Real} <: PointProcess{M}
    pp::P
    tmin::R
    tmax::R
end

"""
    Parameter

An alias for `ComponentVector`, which is how we store point process parameters `θ`.
"""
const Parameter = ComponentVector

Base.eltype(::Type{<:PointProcess{M}}) where {M} = History{M}

StatsBase.params(pp::PointProcess) = ComponentVector(ntfromstruct(pp))

Base.rand(pp::PointProcess, args...) = rand(GLOBAL_RNG, pp, args...)
