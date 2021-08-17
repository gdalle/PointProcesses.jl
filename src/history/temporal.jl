"""
    TemporalHistory{M}

Linear event histories with marks of type `M`.

# Fields

- `times::Vector{Float64}`: vector of event times
- `marks::Vector{M}`: vector of event marks
- `tmin::Float64`: start time
- `tmax::Float64`: end time
"""
mutable struct TemporalHistory{M} <: AbstractHistory{Float64,M}
    times::Vector{Float64}
    marks::Vector{M}
    tmin::Float64
    tmax::Float64
end

function TemporalHistory(times, marks::AbstractMatrix{M}, tmin, tmax) where {M}
    vector_marks = [c for c in eachcol(marks)]
    return TemporalHistory{Vector{M}}(times, vector_marks, tmin, tmax)
end

event_times(h::TemporalHistory) = h.times

event_marks(h::TemporalHistory) = h.marks

mintime(h::TemporalHistory) = h.tmin

maxtime(h::TemporalHistory) = h.tmax

"""
    nb_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)

Count events in `h` during the interval `[tmin, tmax)`.
"""
function nb_events(h::TemporalHistory, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(event_times(h), tmin)
    i_max = searchsortedlast(event_times(h), tmax - eps(tmax))
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
    return maxtime(h) - mintime(h)
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

"""
    append!(h1::TemporalHistory, h2::TemporalHistory)

Add all the events of `h2` at the end of `h1`.
"""
function Base.append!(h1::TemporalHistory, h2::TemporalHistory)
    if has_events(h1) && has_events(h2)
        @assert maximum(event_times(h1)) < minimum(event_times(h2))
    end
    for (t, m) in zip(event_times(h2), event_marks(h2))
        push!(h1, t, m)
    end
    return nothing
end

@doc raw"""
    time_change(h, Λ)

Apply the time rescaling $t \mapsto \Lambda(t)$ to history `h`.
"""
function time_change(h::TemporalHistory, Λ)
    new_times = Λ.(event_times(h))
    new_marks = copy(event_marks(h))
    new_tmin = Λ(mintime(h))
    new_tmax = Λ(maxtime(h))
    return TemporalHistory(new_times, new_marks, new_tmin, new_tmax)
end
