"""
    TemporalHistory{M}

Linear event histories with marks of type `M`.

# Fields

- `times::Vector{Float64}`: vector of event times
- `marks::Vector{M}`: vector of event marks
- `tmin::Float64`: start time
- `tmax::Float64`: end time
"""
@with_kw mutable struct TemporalHistory{M} <: AbstractHistory{Float64,M}
    times::Vector{Float64}
    marks::Vector{M}
    tmin::Float64
    tmax::Float64
end

"""
    nb_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)

Count events in `h` during the interval `[tmin, tmax)`.
"""
function nb_events(h::TemporalHistory, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(h.times, tmin)
    i_max = searchsortedlast(h.times, tmax - eps(tmax))
    return i_max - i_min + 1
end

"""
    has_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)

Check the presence of events in `h` during the interval `[tmin, tmax)`.
"""
function has_events(h::TemporalHistory, tmin = -Inf, tmax = Inf)
    return nb_events(h, tmin, tmax) > 0
end

"""
    duration(h::TemporalHistory)

Compute the difference `h.tmax - h.tmin`.
"""
function duration(h::TemporalHistory)
    return h.tmax - h.tmin
end

"""
    push!(h::TemporalHistory{M}, t::Float64, m::M)

Add event `(t, m)` at the end of history `h`.
"""
function Base.push!(h::TemporalHistory, t, m)
    @assert h.tmin <= t < h.tmax
    push!(h.times, t)
    push!(h.marks, m)
    return nothing
end

@doc raw"""
    time_change(h, Λ)

Apply the time rescaling $t \mapsto \Lambda(t)$ to history `h`.
"""
function time_change(h::TemporalHistory, Λ)
    new_times = Λ.(h.times)
    new_tmin = Λ(h.tmin)
    new_tmax = Λ(h.tmax)
    return TemporalHistory(new_times, h.marks, new_tmin, new_tmax)
end
