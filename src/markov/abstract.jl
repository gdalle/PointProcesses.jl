abstract type AbstractMarkovChain <: AbstractMeasure end

abstract type AbstractMarkovChainPrior <: AbstractMeasure end

function Dists.fit(mctype::Type{<:AbstractMarkovChain}, args...; kwargs...)
    return fit_mle(mctype, args..., kwargs...)
end

function Dists.fit_mle(mctype::Type{<:AbstractMarkovChain}, args...; kwargs...)
    ss = suffstats(mctype, args...; kwargs...)
    return fit_mle(mctype, ss)
end
