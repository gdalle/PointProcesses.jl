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
