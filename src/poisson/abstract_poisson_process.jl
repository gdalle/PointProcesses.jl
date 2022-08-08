abstract type AbstractPoissonProcess{M} <: AbstractPointProcess{M} end

## Defining AbstractPoissonProcess interface

function ground_intensity(pp::AbstractPoissonProcess) end
function mark_distribution(pp::AbstractPoissonProcess) end

function intensity(pp::AbstractPoissonProcess, m)
    return ground_intensity(pp) * densityof(mark_distribution(pp), m)
end

function log_intensity(pp::AbstractPoissonProcess, m)
    return log(ground_intensity(pp)) + logdensityof(mark_distribution(pp), m)
end

## Implementing AbstractPointProcess interface

ground_intensity(pp::AbstractPoissonProcess, t, h) = ground_intensity(pp)
mark_distribution(pp::AbstractPoissonProcess, t, h) = mark_distribution(pp)
intensity(pp::AbstractPoissonProcess, m, t, h) = intensity(pp, m)
log_intensity(pp::AbstractPoissonProcess, m, t, h) = log_intensity(pp, m)

function ground_intensity_bound(pp::AbstractPoissonProcess{M}, t::T, h) where {M,T}
    B = ground_intensity(pp)
    L = typemax(T)
    return (B, L)
end

function integrated_ground_intensity(pp::AbstractPoissonProcess, h, a, b)
    return ground_intensity(pp) * (b - a)
end
