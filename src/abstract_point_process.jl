"""
    AbstractPointProcess{M}

The common supertype of all temporal point processes (i.e. point processes on the real line) with mark type `M`.
"""
abstract type AbstractPointProcess{M} end

@inline DensityInterface.DensityKind(::AbstractPointProcess) = HasDensity()

## Intensity functions

"""
    ground_intensity(pp, h, t)

Compute the ground intensity for a temporal point process `pp` applied to history `h` at time `t`.

The ground intensity quantifies the instantaneous risk of an event with any mark occurring at time `t`:
```
λg(t|h) = Σₘ λ(t,m|h)
```
"""
function ground_intensity(pp::PP, t, h) where {PP<:AbstractPointProcess}
    return error("Not implemented for $PP")
end

"""
    mark_distribution(pp, t, h)

Compute the distribution of marks for a temporal point process `pp` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(pp::PP, t, h) where {PP<:AbstractPointProcess}
    return error("Not implemented for $PP")
end

"""
    intensity(pp, m, t, h)

Compute the conditional intensity for a temporal point process `pp` applied to history `h` and event `(t, m)`.

The conditional intensity function λ(t,m|h) quantifies the instantaneous risk of an event with mark `m` occurring at time `t` after history `h`.
"""
function intensity(pp::AbstractPointProcess, m, t, h)
    return ground_intensity(pp, t, h) * densityof(mark_distribution(pp, t, h), m)
end

"""
    log_intensity(pp, m, t, h)

Compute the logarithm of the conditional intensity for a temporal point process `pp` applied to history `h` and event `(t, m)`.
"""
function log_intensity(pp::AbstractPointProcess, m, t, h)
    return log(ground_intensity(pp, t, h)) + logdensityof(mark_distribution(pp, t, h), m)
end

## Simulation

"""
    ground_intensity_bound(pp, t, h)

Compute a local upper bound on the ground intensity for a temporal point process `pp` applied to history `h` at time `t`.

Return a tuple of the form `(B, L)` satisfying `λg(t|h) ≤ B` for all `u ∈ [t, t+L)`.
"""
function ground_intensity_bound(pp::PP, t, h) where {PP<:AbstractPointProcess}
    return error("Not implemented for $PP")
end

## Learning

"""
    integrated_ground_intensity(pp, h, a, b)

Compute the integrated ground intensity (or compensator) `Λ(t|h)` for a temporal point process `pp` applied to history `h` on interval `[a, b]`:
```
Λ(h) = ∫ λg(t|h) dt
```
"""
function integrated_ground_intensity(pp::PP, h, a, b) where {PP<:AbstractPointProcess}
    return error("Not implemented for $PP")
end

"""
    logdensityof(pp, h)

Compute the log probability density function for a temporal point process `pp` applied to history `h`:
```
ℓ(h) = Σₖ log λ(tₖ|hₖ) - Λ(h)
```
The default method uses a loop over events combined with [`integrated_ground_intensity`](@ref), but it should be reimplemented for specific processes if faster computation is possible.
"""
function DensityInterface.logdensityof(pp::AbstractPointProcess, h)
    l = -integrated_ground_intensity(pp, h, min_time(h), max_time(h))
    for (t, m) in zip(event_times(h), event_marks(h))
        l += log_intensity(pp, m, t, h)
    end
    return l
end

"""
    check_residuals(pp, h)

Check whether the point process `pp` is a good fit for history `h` by applying Ogata's time rescaling method.
If `(tₖ)ₖ` is a temporal point process with intensity `λ`, then `(Λ(tₖ))ₖ` is a standard temporal Poisson process.
"""
function check_residuals(pp::AbstractPointProcess, h)
    return error("Work in progress.")
end
