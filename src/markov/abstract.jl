abstract type AbstractMarkovChain end

abstract type AbstractMarkovChainPrior end

Base.rand(mc::AbstractMarkovChain, args...) = rand(Random.GLOBAL_RNG, mc, args...)

function Distributions.fit(mctype::Type{<:AbstractMarkovChain}, args...)
    return fit_mle(mctype, args...)
end

function Distributions.fit_mle(mctype::Type{<:AbstractMarkovChain}, args...)
    ss = suffstats(mctype, args...)
    return fit_mle(mctype, ss)
end
