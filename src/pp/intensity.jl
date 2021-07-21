
@doc raw"""
    intensity(pptype, θ, h, t, m)

Compute the conditional intensity for a point process of type `pptype` with parameter `θ` applied to history `h` and event `(t, m)`.

The conditional intensity function quantifies the instantaneous risk  of an event with mark `m` occurring at time `t`[^Rasmussen_2018].

[^Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].
"""
function intensity(pptype, θ, h, t, m) end

"""
    intensity(pp, h, t, m)

Fall back on [`intensity(pptype, θ, h, t, m)`](@ref).
"""
function intensity(pp::PointProcess{M}, h::History{M}, t::Real, m::M) where {M}
    return intensity(typeof(pp), params(pp), h, t, m)
end

@doc raw"""
    mark_distribution(pptype, θ, h, t)

Compute the distribution of marks for a point process of type `pptype` with parameter `θ` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(pptype, θ, h, t) end

"""
    mark_distribution(pp, h, t)

Fall back on [`mark_distribution(pptype, θ, h, t)`](@ref).
"""
function mark_distribution(pp::PointProcess{M}, h::History{M}, t::Real) where {M}
    return mark_distribution(typeof(pp), params(pp), h, t)
end

@doc raw"""
    ground_intensity(pptype, θ, h, t)

Compute the ground intensity for a point process of type `pptype` with parameter `θ` applied to history `h` at time `t`.

The ground intensity quantifies the instantaneous risk of an event with any mark occurring at time `t`[^Rasmussen_2018]. It can be expressed as
```math
    \lambda_g(t|h) = \sum_{m \in \mathcal{M}} \lambda(t, m|h).
```
"""
function ground_intensity(pptype, θ, h, t) end

"""
    ground_intensity(pp, h, t)

Fall back on [`ground_intensity(pptype, θ, h, t)`](@ref).
"""
function ground_intensity(pp::PointProcess{M}, h::History{M}, t::Real) where {M}
    return ground_intensity(typeof(pp), params(pp), h, t)
end

@doc raw"""
    ground_intensity_bound(pptype, θ, h, t)

Compute a local upper bound on the ground intensity for a point process of type `pptype` with parameter `θ` applied to history `h` at time `t`[^Rasmussen_2018].

Return a tuple of the form $(B, L)$ satisfying
```math
    \forall u \in [t, t+L), \quad \lambda_g(t|h) \leq B
```
"""
function ground_intensity_bound(pptype, θ, h, t) end

"""
    ground_intensity_bound(pp, h, t)

Fall back on [`ground_intensity_bound(pptype, θ, h, t)`](@ref).
"""
function ground_intensity_bound(pp::PointProcess{M}, h::History{M}, t::Real) where {M}
    return ground_intensity_bound(typeof(pp), params(pp), h, t)
end
