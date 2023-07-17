Base.@kwdef struct MarkedPoissonProcessStats{R1<:Real,R2<:Real,DSS}
    nb_events::R1
    duration::R2
    mark_dist_suffstats::DSS
end

## Compute sufficient stats

function Distributions.suffstats(
    ::Type{MarkedPoissonProcess{R,M,D}}, h::History
) where {R,M,D}
    return MarkedPoissonProcessStats(;
        nb_events=nb_events(h),
        duration=duration(h),
        mark_dist_suffstats=suffstats(D, event_marks(h)),
    )
end

function Distributions.suffstats(
    ::Type{MarkedPoissonProcess{R,M,D}},
    histories::AbstractVector{<:History},
    weights::AbstractVector{W},
) where {R,M,D,W}
    nb_events_total = sum(nb_events, histories)
    duration_total = sum(duration, histories)
    event_marks_concat = mapreduce(event_marks, vcat, histories)
    weights_concat = reduce(
        vcat, (fill(w, nb_events(h)) for (w, h) in zip(weights, histories))
    )
    return MarkedPoissonProcessStats(;
        nb_events=nb_events_total,
        duration=duration_total,
        mark_dist_suffstats=suffstats(D, event_marks_concat, weights_concat),
    )
end

function Distributions.suffstats(
    ::Type{MarkedPoissonProcess{R,M,D}}, histories::AbstractVector{<:History}
) where {R,M,D}
    return suffstats(MarkedPoissonProcess{R,M,D}, histories, ones(length(histories)))
end
