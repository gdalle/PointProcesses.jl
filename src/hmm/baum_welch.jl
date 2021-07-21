function expectation_maximization_step!(
    logα::Matrix,
    logβ::Matrix,
    logγ::Matrix,
    logξ::Array,
    obs_logpdf::Matrix,
    hmm::HMM,
    observations::Vector,
)
    update_observation_likelihood!(obs_logpdf, hmm, observations)
    logL = forward_backward_log!(logα, logβ, logγ, logξ, obs_logpdf, hmm)
    hmm.transitions = fit(typeof(hmm.transitions), logγ, logξ)
    hmm.emissions = [
        fit(typeof(hmm.emissions[s]), @view logγ[:, s], observations) for
        s = 1:nstates(hmm)
    ]
end

function expectation_maximization(hmm::HMM, observations::Vector)

end
