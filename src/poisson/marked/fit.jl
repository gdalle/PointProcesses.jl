## Fit from sufficient stats

function StatsAPI.fit(
    ::Type{MarkedPoissonProcess{R,M,D}}, ss::MarkedPoissonProcessStats
) where {R,M,D}
    λ = convert(R, ss.nb_events / ss.duration)
    mark_dist = fit(D, ss.mark_dist_suffstats)
    return MarkedPoissonProcess(λ, mark_dist)
end

## Fit from observations

function StatsAPI.fit(
    pptype::Type{MarkedPoissonProcess{R,M,D}}, args...; kwargs...
) where {R,M,D}
    ss = suffstats(pptype, args...; kwargs...)
    return fit(pptype, ss)
end
