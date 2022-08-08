"""
    BoundedPointProcess{M,P,T}

Store a temporal point process `P` with pre-defined start and end times.
"""
struct BoundedPointProcess{M,P<:AbstractPointProcess{M},T<:Real} <: AbstractPointProcess{M}
    pp::P
    tmin::T
    tmax::T
end

min_time(pp::BoundedPointProcess) = pp.tmin
max_time(pp::BoundedPointProcess) = pp.tmax

function Base.rand(rng::AbstractRNG, pp::BoundedPointProcess)
    return rand(rng, pp.pp, min_time(pp), max_time(pp))
end
