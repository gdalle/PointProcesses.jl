struct MultivariatePoissonProcessStats{R1<:Real,R2<:Real}
    nb_events::Vector{R1}
    duration::R2
end

function add_suffstats(
    ss1::MultivariatePoissonProcessStats{R1,R2}, ss2::MultivariatePoissonProcessStats{R1,R2}
) where {R1<:Real,R2<:Real}
    nb_events = ss1.nb_events .+ ss2.nb_events
    duration = ss1.duration + ss2.duration
    return MultivariatePoissonProcessStats{R1,R2}(nb_events, duration)
end

## Compute sufficient stats

function Distributions.suffstats(
    ::Type{MultivariatePoissonProcess{R}}, history::History{<:Integer,<:Real}
) where {R}
    M = maximum_mark(history)
    nb_events = zeros(Int, M)
    for m in event_marks(history)
        nb_events[m] += 1
    end
    return MultivariatePoissonProcessStats(nb_events, duration(history))
end

function Distributions.suffstats(
    ::Type{MultivariatePoissonProcess{R}},
    histories::AbstractVector{<:History{<:Integer,<:Real}},
    weights::AbstractVector{W},
) where {R,W<:Real}
    M = mapreduce(maximum_mark, max, histories)
    total_nb_events = zeros(W, M)
    total_duration = zero(W)
    for (history, weight) in zip(histories, weights)
        total_duration += weight * duration(history)
        for m in event_marks(history)
            total_nb_events[m] += weight
        end
    end
    return MultivariatePoissonProcessStats(total_nb_events, total_duration)
end
