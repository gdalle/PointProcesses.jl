function forward_nolog!(α::Matrix, c::Vector, obs_pdf::Matrix, hmm::HiddenMarkovModel)
    T, S = length(obs_logpdf), nstates(hmm)
    π0, P = hmm.transitions.π0, hmm.transitions.P
    for i = 1:S
        α[1, i] = π0[i] * obs_pdf[1, i]
    end
    c[1] = 1 / sum(@view α[1, :])  # scaling
    for i = 1:S
        α[1, i] *= c[1]
    end
    for t = 1:T-1
        for j = 1:S
            α[t+1, j] = sum(α[t, i] * P[i, j]) * obs_pdf[t+1, j]
        end
        c[t+1] = 1 / sum(@view α[t+1, :])  # scaling
        for j = 1:S
            α[t+1, j] *= c[t+1]
        end
    end
    for t = 1:T
        if all_zeros(@view α[t, :])
            throw(OverflowError("Probabilities are too small in forward step."))
        end
    end
    return logα
end

function forward_log!(logα::Matrix, obs_logpdf::Matrix, hmm::HiddenMarkovModel)
    T, S = length(obs_logpdf), nstates(hmm)
    logπ0, logP = log.(hmm.transitions.π0), log.(hmm.transitions.P)
    for i = 1:S
        logα[1, i] = logπ0[i] + obs_logpdf[1, i]
    end
    for t = 1:T-1
        for j = 1:S
            logα[t+1, j] = logsumexp(logα[t, :] + logP[:, s]) + obs_logpdf[t+1, j]  # TODO
        end
    end
    for t = 1:T
        if all_nan(@view logα[t, :])
            throw(OverflowError("Log probabilities are too small in forward step."))
        end
    end
    return logα
end

function backward_nolog!(β::Matrix, c::Vector, obs_pdf::Matrix, hmm::HiddenMarkovModel)
    T, S = length(obs_logpdf), nstates(hmm)
    π0, P = hmm.transitions.π0, hmm.transitions.P

    for i = 1:S
        β[T, i] = 0.
    end
    for t = T-1:-1:1
        for s = 1:S
            β[t, i] = sum(P[i, j] * obs_pdf[t+1, j] * β[t+1, j]) * c[t]
        end
    end
    for t = 1:T
        if all_zeros(@view β[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
end


function backward_log!(logβ::Matrix, obs_logpdf::Matrix, hmm::HiddenMarkovModel)
    T, S = length(obs_logpdf), nstates(hmm)
    logπ0, logP = log.(hmm.transitions.π0), log.(hmm.transitions.P)

    for i = 1:S
        logβ[T, i] = 0.
    end
    for t = T-1:-1:1
        for i = 1:i
            logβ[t, i] = logsumexp(logP[s, :] .+ obs_logpdf[t+1, :] .+ logβ[t+1, :])  # TODO
        end
    end
    for t = 1:T
        if all_nan(@view logβ[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
end

function forward_backward_nolog!(
    α::Matrix,
    c::Vector,
    β::Matrix,
    γ::Matrix,
    ξ::Array,
    obs_pdf::Matrix,
    hmm::HiddenMarkovModel,
)
    T, S = length(observations), nstates(hmm.transitions)

    forward_nolog!(α, c, hmm, obs_pdf)
    backward_nolog!(β, c, hmm, obs_pdf)

    for t = 1:T
        for i = 1:S
            γ[t, i] = α[t, i] * β[t, i]
        end
        sumγ = sum(@view γ[t, :])
        for i = 1:S
            γ[t, i] /= sumγ
        end
    end

    for t = 1:T-1
        for i = 1:S, j = 1:S
            ξ[t, i, j] = α[t, i] * P[i, j] * obs_pdf[t+1, j] * β[t+1, j]
        end
        sumξ = sum(@view ξ[t, :, :])
        for i = 1:S, j = 1:S
            ξ[t, i, j] /= sumξ
        end
    end

    L = sum(@view α[T, :])

    return L
end


function forward_backward_log!(
    logα::Matrix,
    logβ::Matrix,
    logγ::Matrix,
    logξ::Array,
    obs_logpdf::Matrix,
    hmm::HiddenMarkovModel,
)
    T, S = length(observations), nstates(hmm.transitions)
    logπ0, logP = log.(hmm.transitions.π0), log.(hmm.transitions.P)

    forward_log!(α, obs_logpdf, hmm)
    backward_log!(β, obs_logpdf, hmm)

    for t = 1:T
        for i = 1:S
            logγ[t, i] = logα[t, i] + logβ[t, i] - lse
        end
        logsumγ = logsumexp(@view logγ[t, :])
        for i = 1:S
            logγ[t, i] -= logsumγ
        end
    end

    for t = 1:T-1
        for i = 1:S, j = 1:S
            logξ[t, i, j] = logα[t, i] + logP[i, j] + obs_logpdf[t+1, j] + logβ[t+1, j]
        end
        logsumξ = logsumexp(@view logξ[t, :, :])
        for i = 1:S, j = 1:S
            logξ[t, i, j] -= logsumξ
        end
    end

    logL = logsumexp(logα[T, :])

    return logL
end
