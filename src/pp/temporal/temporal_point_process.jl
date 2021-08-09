"""
    TemporalPointProcess{M}

The common supertype of all point processes on the real line with mark type `M`.
"""
abstract type TemporalPointProcess{M} <: AbstractPointProcess{Float64,M} end

"""
    Bounded{M, P}

Store a temporal point process `P` with pre-defined start and end times.
"""
struct Bounded{M, P<:TemporalPointProcess{M}} <: TemporalPointProcess{M}
    pp::P
    tmin::Float64
    tmax::Float64
end
