## Fit from sufficient stats

function Distributions.fit_mle(
    ::Type{MarkedPoissonProcess{R,M,D}}, ss::MarkedPoissonProcessStats
) where {R,M,D}
    λ = convert(R, ss.nb_events / ss.duration)
    mark_dist = fit_mle(D, ss.mark_dist_suffstats)
    return MarkedPoissonProcess(λ, mark_dist)
end

## Fit from observations

function Distributions.fit_mle(
    pptype::Type{MarkedPoissonProcess{R,M,D}}, args...; kwargs...
) where {R,M,D}
    ss = suffstats(pptype, args...; kwargs...)
    return fit_mle(pptype, ss)
end
