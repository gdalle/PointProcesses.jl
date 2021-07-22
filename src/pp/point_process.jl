"""
    PointProcess{M}

The common supertype of all point processes with mark type `M`.
"""
abstract type PointProcess{M} end

"""
    TimedPointProcess{M,P<:PointProcess{M}}

Abstract point process with built-in start and end times.

# Fields
- `pp::P`: underlying point process
- `tmin::Float64`: start time
- `tmax::Float64`: end time
"""
struct TimedPointProcess{M,P<:PointProcess{M}} <: PointProcess{M}
    pp::P
    tmin::Float64
    tmax::Float64
end

"""
    Parameter

An alias for `ComponentVector`, which is how we store point process parameters `θ`.
"""
const Parameter = ComponentVector

StatsBase.params(pp::PointProcess) = ComponentVector(ntfromstruct(pp))
