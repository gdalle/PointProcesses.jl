abstract type AbstractMarkovChain end

StatsBase.params(mc::AbstractMarkovChain) = ntfromstruct(mc)

nstates(mc::AbstractMarkovChain) = length(mc.π0)

function Distributions.fit(mctype::Type{<:AbstractMarkovChain}, args...)
    return fit_mle(mctype, args...)
end

function Distributions.fit_mle(mctype::Type{<:AbstractMarkovChain}, args...)
    ss = suffstats(mctype, args...)
    return fit_mle(mctype, ss)
end
