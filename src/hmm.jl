struct HMM{E}
    π0::Vector{Float64}
    transitions::DiscreteMarkovChain
    emissions::E
end

nstates(hmm::HMM) = length(hmm.emissions)

# Simulation

function Distributions.rand(hmm::HMM, T)
    states = rand(hmm.transitions, T)
    observations = [rand(hmm.emissions[states[t]]) for t = 1:T]
    return states, observations
end
