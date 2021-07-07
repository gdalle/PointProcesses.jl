import Base: push!, append!

mutable struct History{MarkType}
    times
    marks::Vector{MarkType}
    tmin
    tmax
end

function push!(history::History, t, m)
    @assert t >= history.tmax
    push!(history.times, t)
    push!(history.marks, m)
    history.tmax = max(t, history.tmax)
end

function append!(history1::History, history2::History)
    @assert history1.tmax == history2.tmin
    append!(history1.times, history2.times)
    append!(history1.marks, history2.marks)
    history1.tmax = history2.tmax
end

function nb_eventimes(history, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(history.times, tmin)
    i_max = searchsortedlast(history.times, tmax - eps())
    return i_max - i_min + 1
end

function has_eventimes(history, tmin = -Inf, tmax = Inf)
    return nb_eventimes(history, tmin, tmax) > 0
end

