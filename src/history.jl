mutable struct History{M}
    times::Vector{Float64}
    marks::Vector{M}
    tmin::Float64
    tmax::Float64
end

function nb_events(history, tmin = -Inf, tmax = Inf)
    i_min = searchsortedfirst(history.times, tmin)
    i_max = searchsortedlast(history.times, tmax - eps())
    return i_max - i_min + 1
end

function has_events(history, tmin = -Inf, tmax = Inf)
    return nb_events(history, tmin, tmax) > 0
end

function duration(history)
    return history.tmax - history.tmin
end

function Base.push!(history::History{M}, t::Float64, m::M) where {M}
    push!(history.times, t)
    push!(history.marks, m)
end
