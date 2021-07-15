import Base: push!, append!

mutable struct History{M}
    times::Vector{Float64}
    marks::Vector{M}
    tmin::Float64
    tmax::Float64
end

function get_times(history)
    return history.times
end

function get_marks(history)
    return history.marks
end

function get_tmin(history)
    return history.tmin
end

function get_tmax(history)
    return history.tmax
end

function set_tmin!(history, tmin)
    history.tmin = tmin
end

function set_tmax!(history, tmax)
    history.tmax = tmax
end

function push!(history::History{M}, t::Float64, m::M) where {M}
    push!(history.times, t)
    push!(history.marks, m)
end

function append!(history1::History{M}, history2::History{M}) where {M}
    @assert history1.tmax == history2.tmin
    append!(history1.times, history2.times)
    append!(history1.marks, history2.marks)
    set_tmax!(history1, history2.tmax)
end

function nb_events(history, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(history.times, tmin)
    i_max = searchsortedlast(history.times, tmax - eps())
    return i_max - i_min + 1
end

function has_events(history, tmin = -Inf, tmax = Inf)
    return nb_events(history, tmin, tmax) > 0
end

