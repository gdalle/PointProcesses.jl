"""
    BoundedPointProcess{M,P,T}

Store a temporal point process `P` with pre-defined start and end times.
"""
struct BoundedPointProcess{M,P<:AbstractPointProcess{M},T<:Real} <: AbstractPointProcess{M}
    pp::P
    tmin::T
    tmax::T
end

min_time(bpp::BoundedPointProcess) = bpp.tmin
max_time(bpp::BoundedPointProcess) = bpp.tmax

function Base.rand(rng::AbstractRNG, bpp::BoundedPointProcess)
    return rand(rng, bpp.pp, min_time(bpp), max_time(bpp))
end

ground_intensity(bpp::BoundedPointProcess, t, h) = ground_intensity(bpp.pp, t, h)
mark_distribution(bpp::BoundedPointProcess, t, h) = mark_distribution(bpp.pp, t, h)
intensity(bpp::BoundedPointProcess, m, t, h) = intensity(bpp.pp, m, t, h)
log_intensity(bpp::BoundedPointProcess, m, t, h) = log_intensity(bpp.pp, m, t, h)

function ground_intensity_bound(bpp::BoundedPointProcess, t, h)
    return ground_intensity_bound(bpp.pp, t, h)
end

function integrated_ground_intensity(bpp::BoundedPointProcess, h, a, b)
    return integrated_ground_intensity(bpp.pp, h, a, b)
end
