"""
    TemporalPointProcess{L<:Real, M}

The common supertype of all temporal point processes (i.e. point processes on the real line) with mark type `M`.
"""
abstract type TemporalPointProcess{M} <: AbstractPointProcess{Float64,M} end

"""
    BoundedTemporalPointProcess{L<:Real, M}

Store a temporal point process `P` with pre-defined start and end times.
"""
@with_kw struct BoundedTemporalPointProcess{M,P<:TemporalPointProcess{M}} <: TemporalPointProcess{M}
    pp::P
    tmin::Float64
    tmax::Float64
end
