"""
    TemporalPointProcess{M,T}

The common supertype of all temporal point processes (i.e. point processes on the real line) with mark type `M` and time type `T`.
"""
abstract type TemporalPointProcess{M,T<:Real} end

@inline DensityInterface.DensityKind(::TemporalPointProcess) = HasDensity()

"""
    BoundedTemporalPointProcess{M,T}

Store a temporal point process `P` with pre-defined start and end times.
"""
struct BoundedTemporalPointProcess{M,T<:Real,P<:TemporalPointProcess{M,T}} <:
       TemporalPointProcess{M,T}
    pp::P
    tmin::T
    tmax::T
end

## Intensity functions

@doc raw"""
    intensity(pp, m, t, h)

Compute the conditional intensity for a temporal point process `pp` applied to history `h` and event `(t, m)`.

The conditional intensity function $\lambda(t, m | h)$ quantifies the instantaneous risk  of an event with mark `m` occurring at time `t`.
"""
function intensity(pp::TemporalPointProcess, m, t, h) end

@doc raw"""
    log_intensity(pp, m, t, h)

Compute the logarithm of the conditional intensity for a temporal point process `pp` applied to history `h` and event `(t, m)`.
"""
log_intensity(pp::TemporalPointProcess, m, t, h) = log(intensity(pp, m, t, h))

@doc raw"""
    mark_distribution(pp, t, h)

Compute the distribution of marks for a temporal point process `pp` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(pp::TemporalPointProcess, t, h) end

@doc raw"""
    ground_intensity(pp, h, t)

Compute the ground intensity for a temporal point process `pp` applied to history `h` at time `t`.

The ground intensity quantifies the instantaneous risk of an event with any mark occurring at time `t`. It can be expressed as
```math
    \lambda_g(t|h) = \sum_{m \in \mathcal} \lambda(t, m|h).
```
"""
function ground_intensity(pp::TemporalPointProcess, t, h) end

@doc raw"""
    ground_intensity_bound(pp, t, h)

Compute a local upper bound on the ground intensity for a temporal point process `pp` applied to history `h` at time `t`.

Return a tuple of the form $(B, L)$ satisfying
```math
    \forall u \in [t, t+L), \quad \lambda_g(t|h) \leq B
```
"""
function ground_intensity_bound(pp::TemporalPointProcess, t, h) end

## Learning

@doc raw"""
    logdensityof(pp, h)

Compute the log probability density function for a temporal point process `pp` applied to history `h`:
```math
    \log f(h) = \sum_i \log \lambda(t_i | h_i) - \Lambda(h)
```

The default method uses a loop over events combined with [`integrated_ground_intensity`](@ref), but it should be reimplemented for specific processes if faster computation is possible.
"""
function DensityInterface.logdensityof(pp::TemporalPointProcess, h::History)
    l = -integrated_ground_intensity(pp, h)
    for (t, m) in zip(event_times(h), event_marks(h))
        l += log_intensity(pp, m, t, h)
    end
    return l
end
