@doc raw"""
    intensity(pp, h, t, m)

Compute the conditional intensity for a temporal point process `pp` applied to history `h` and event `(t, m)`.

The conditional intensity function ``\lambda(t, m | h)`` quantifies the instantaneous risk  of an event with mark `m` occurring at time `t`[^Rasmussen_2018].

[^Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].
"""
function intensity(pp::TemporalPointProcess, h::TemporalHistory, t, m)
    error("not implemented")
end

@doc raw"""
    mark_distribution(pp, h, t)

Compute the distribution of marks for a temporal point process `pp` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(pp::TemporalPointProcess, h::TemporalHistory, t)
    error("not implemented")
end

@doc raw"""
    ground_intensity(pp, h, t)

Compute the ground intensity for a temporal point process `pp` applied to history `h` at time `t`.

The ground intensity quantifies the instantaneous risk of an event with any mark occurring at time `t`[^Rasmussen_2018]. It can be expressed as
```math
    \lambda_g(t|h) = \sum_{m \in \mathcal} \lambda(t, m|h).
```
"""
function ground_intensity(pp::TemporalPointProcess, h::TemporalHistory, t)
    error("not implemented")
end

@doc raw"""
    ground_intensity_bound(pp, θ, h, t)

Compute a local upper bound on the ground intensity for a temporal point process `pp` applied to history `h` at time `t`[^Rasmussen_2018].

Return a tuple of the form $(B, L)$ satisfying
```math
    \forall u \in [t, t+L), \quad \lambda_g(t|h) \leq B
```
"""
function ground_intensity_bound(pp::TemporalPointProcess, h::TemporalHistory, t)
    error("not implemented")
end
