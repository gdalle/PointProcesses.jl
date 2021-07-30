"""
    History{M, R}

A container for linear event histories with times of type `R` and marks of type `M`.

# Fields

- `times::Vector{R}`: vector of event times
- `marks::Vector{M}`: vector of event marks
- `tmin::R`: start time
- `tmax::R`: end time

# Examples

```jldoctest
julia> h = History([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0)
History{String}([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0)

julia> duration(h)
2.0

julia> nb_events(h)
3

julia> nb_events(h, 1.0, 2.0)
1

julia> has_events(h)
true

julia> has_events(h, 1.5, 2.0)
false

julia> push!(h, 1.7, "d")

julia> has_events(h, 1.5, 2.0)
true
```
"""
mutable struct History{M, R}
    times::Vector{R}
    marks::Vector{M}
    tmin::R
    tmax::R
end

"""
    nb_events(h::History, tmin=-Inf, tmax=Inf)

Count events in `h` during the interval `[tmin, tmax)`.
"""
function nb_events(h::History, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(h.times, tmin)
    i_max = searchsortedlast(h.times, tmax - eps(tmax))
    return i_max - i_min + 1
end

"""
    has_events(h::History, tmin=-Inf, tmax=Inf)

Check the presence of events in `h` during the interval `[tmin, tmax)`.
"""
function has_events(h::History, tmin = -Inf, tmax = Inf)
    return nb_events(h, tmin, tmax) > 0
end

"""
    duration(h::History)

Compute the difference `h.tmax - h.tmin`.
"""
function duration(h::History)
    return h.tmax - h.tmin
end

"""
    push!(h::History{M}, t::Float64, m::M)

Add event `(t, m)` at the end of history `h`.
"""
function Base.push!(h::History, t, m)
    @assert h.tmin <= t < h.tmax
    push!(h.times, t)
    push!(h.marks, m)
    return nothing
end

@doc raw"""
    time_change(h, Λ)

Apply the time rescaling $t \mapsto \Lambda(t)$ to history `h`.
"""
function time_change(h::History, Λ)
    new_times = Λ.(h.times)
    new_tmin = Λ(h.tmin)
    new_tmax = Λ(h.tmax)
    return History(new_times, h.marks, new_tmin, new_tmax)
end
