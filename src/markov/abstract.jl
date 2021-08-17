abstract type AbstractMarkovChain end

abstract type AbstractMarkovChainPrior end

function Distributions.fit(mctype::Type{<:AbstractMarkovChain}, args...; kwargs...)
    return fit_mle(mctype, args..., kwargs...)
end

function Distributions.fit_mle(mctype::Type{<:AbstractMarkovChain}, args...; kwargs...)
    ss = suffstats(mctype, args...; kwargs...)
    return fit_mle(mctype, ss)
end
