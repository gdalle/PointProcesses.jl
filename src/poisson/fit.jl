## Fit MLE

function StatsAPI.fit(
    ::Type{PoissonProcess{R,D}}, ss::PoissonProcessStats{R1,R2}
) where {R<:Real,D,R1<:Real,R2<:Real}
    λ = convert(R, ss.nb_events / ss.duration)
    mark_dist = fit(D, ss.marks, ss.weights)
    return PoissonProcess(λ, mark_dist)
end

function StatsAPI.fit(pptype::Type{<:PoissonProcess}, args...; kwargs...)
    ss = suffstats(pptype, args...; kwargs...)
    return fit(pptype, ss)
end

## Bayesian fit (only for MultivariatePoissonProcess)
# TODO: Replace PoissonProcess{R,Categorical{R,Vector{R}}} by
# MultivariatePoissonProcess{R}
# when everything else is OK

function fit_map(
    pptype::Type{PoissonProcess{R,Categorical{R,Vector{R}}}},
    prior::PoissonProcessPrior,
    ss::PoissonProcessStats,
) where {R<:Real}
    (; α, β) = prior
    posterior_nb_events = [sum(ss.marks .== i) for i in 1:length(α)] .+ α
    posterior_duration = ss.duration + β
    λ = convert(Vector{R}, posterior_nb_events ./ posterior_duration)
    return pptype(sum(λ), Categorical(λ ./ sum(λ)))
end

function fit_map(
    pptype::Type{<:PoissonProcess{R,Categorical{R,Vector{R}}} where {R<:Real}},
    prior::PoissonProcessPrior,
    args...;
    kwargs...,
)
    ss = suffstats(pptype, args..., kwargs...)
    return fit_map(pptype, prior, ss)
end
