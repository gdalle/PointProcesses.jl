function StatsAPI.fit(
    ::Type{MarkedPoissonProcess{M,R,D}},
    histories::AbstractVector{<:History},
    weights::AbstractVector{<:Real};
) where {R,M,D}
    nb_events_total = sum(nb_events, histories)
    duration_total = sum(duration, histories)
    event_marks_concat = mapreduce(event_marks, vcat, histories)
    weights_concat = reduce(
        vcat, (fill(w, nb_events(h)) for (w, h) in zip(weights, histories))
    )
    λ = convert(R, nb_events_total / duration_total)
    mark_dist = fit(D, event_marks_concat, weights_concat)
    return MarkedPoissonProcess{M}(λ, mark_dist)
end

function StatsAPI.fit(
    pptype::Type{<:MarkedPoissonProcess}, histories::AbstractVector{<:History}
)
    weights = ones(length(histories))
    return fit(pptype, histories, weights)
end

function StatsAPI.fit(pptype::Type{<:MarkedPoissonProcess}, h::History)
    return fit(pptype, [h])
end
