## Fit from sufficient stats

function StatsAPI.fit(
    ::Type{PoissonProcess{M,R,D}}, ss::PoissonProcessStats{R1,R2,M}
) where {R,M,D,R1,R2}
    λ = convert(R, ss.nb_events ./ ss.duration)
    mark_dist = fit(D, ss.marks, ss.weights)
    return PoissonProcess(λ, mark_dist)
end

function fit_map(
    pptype::Type{PoissonProcess{M,R,D}},
    prior::PoissonProcessPrior,
    ss::PoissonProcessStats{R1,R2,M},
) where {R,M,D,R1,R2}
    (; α, β) = prior
    posterior_nb_events = ss.nb_events .+ α 
    ## Note that there was a -1 here but I think it was a mistake
    posterior_duration = ss.duration + β
    ss_posterior = PoissonProcessStats(
        posterior_nb_events, posterior_duration, ss.marks, ss.weights
    )
    return fit(pptype, ss_posterior)
end

## Fit from observations

function StatsAPI.fit(
    pptype::Type{PoissonProcess{M,R,D}}, args...; kwargs...
) where {M,R,D}
    ss = suffstats(pptype, args...; kwargs...)
    return fit(pptype, ss)
end

function fit_map(
    pptype::Type{PoissonProcess{M,R,D}},
    prior::PoissonProcessPrior,
    args...;
    kwargs...,
) where {M,R,D}
    ss = suffstats(pptype, args..., kwargs...)
    return fit_map(pptype, prior, ss)
end
