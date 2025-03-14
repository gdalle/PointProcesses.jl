"""
    AbstractPoissonProcess <: AbstractPointProcess

Common interface for all temporal Poisson processes, that is, temporal point processes for which the intensity is not a function of past history.

Implements some fallbacks for the `AbstractPointProcess` interface which accept fewer arguments.
"""
abstract type AbstractPoissonProcess <: AbstractPointProcess end

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
mark_distribution(pp::AbstractPoissonProcess, t) = mark_distribution(pp) # For simulate_ogata
intensity(pp::AbstractPoissonProcess, m, t, h) = intensity(pp, m)
log_intensity(pp::AbstractPoissonProcess, m, t, h) = log_intensity(pp, m)

function ground_intensity_bound(pp::AbstractPoissonProcess, t::T, h) where {T}
    B = ground_intensity(pp)
    L = typemax(T)
    return (B, L)
end

function integrated_ground_intensity(pp::AbstractPoissonProcess, h, a, b)
    return ground_intensity(pp) * (b - a)
end
