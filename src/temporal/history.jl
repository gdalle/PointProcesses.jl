"""
    History{M,T<:Real}

Linear event histories with marks of type `M` and locations of real type `T`.

# Fields

- `times::Vector{T}`: vector of event times
- `marks::Vector{M}`: vector of event marks
- `tmin::T`: start time
- `tmax::T`: end time
"""
Base.@kwdef mutable struct History{M,T<:Real}
    times::Vector{T}
    marks::Vector{M}
    tmin::T
    tmax::T
end

event_times(h::History) = h.times

event_marks(h::History) = h.marks

"""
    min_time(h)

Return the starting time of `h` (not the same as the first event time).
"""
min_time(h::History) = h.tmin

"""
    max_time(h)

Return the end time of `h` (not the same as the last event time).
"""
max_time(h::History) = h.tmax

"""
    nb_events(h::History, tmin=-Inf, tmax=Inf)

Count events in `h` during the interval `[tmin, tmax)`.
"""
function nb_events(h::History, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(event_times(h), tmin)
    i_max = searchsortedlast(event_times(h), tmax - eps(tmax))
    return i_max - i_min + 1
end

"""
    has_events(h::History, tmin=-Inf, tmax=Inf)

Check the presence of events in `h` during the interval `[tmin, tmax)`.
"""
has_events(h::History, tmin = -Inf, tmax = Inf) = nb_events(h, tmin, tmax) > 0

"""
    duration(h::History)

Compute the difference `h.tmax - h.tmin`.
"""
duration(h::History) = max_time(h) - min_time(h)

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

"""
    append!(h1::History, h2::History)

Add all the events of `h2` at the end of `h1`.
"""
function Base.append!(h1::History, h2::History)
    max_time(h1) ≈ min_time(h2) || return false
    append!(h1.times, h2.times)
    append!(h1.marks, h2.marks)
    h1.tmax = h2.tmax
    return true
end

@doc raw"""
    time_change(h, Λ)

Apply the time rescaling $t \mapsto \Lambda(t)$ to history `h`.
"""
function time_change(h::History, Λ)
    new_times = Λ.(event_times(h))
    new_marks = copy(event_marks(h))
    new_tmin = Λ(min_time(h))
    new_tmax = Λ(max_time(h))
    return History(times=new_times, marks=new_marks, tmin=new_tmin, tmax=new_tmax)
end
