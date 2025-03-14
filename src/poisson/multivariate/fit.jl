## Fit from sufficient stats

function StatsAPI.fit(
    ::Type{MultivariatePoissonProcess{R}}, ss::MultivariatePoissonProcessStats
) where {R}
    λ = convert(Vector{R}, ss.nb_events ./ ss.duration)
    return MultivariatePoissonProcess{R}(λ)
end

function fit_map(
    pptype::Type{MultivariatePoissonProcess{R}},
    prior::MultivariatePoissonProcessPrior,
    ss::MultivariatePoissonProcessStats,
) where {R}
    (; λ_α, λ_β) = prior
    posterior_nb_events = ss.nb_events .+ λ_α
    ## Note that there was a -1 above but I think it was a mistake
    posterior_duration = ss.duration + λ_β
    ss_posterior = MultivariatePoissonProcessStats(posterior_nb_events, posterior_duration)
    return fit(pptype, ss_posterior)
end

## Fit from observations

function StatsAPI.fit(
    pptype::Type{MultivariatePoissonProcess{R}}, args...; kwargs...
) where {R}
    ss = suffstats(pptype, args...; kwargs...)
    return fit(pptype, ss)
end

function fit_map(
    pptype::Type{MultivariatePoissonProcess{R}},
    prior::MultivariatePoissonProcessPrior,
    args...;
    kwargs...,
) where {R}
    ss = suffstats(pptype, args..., kwargs...)
    return fit_map(pptype, prior, ss)
end
