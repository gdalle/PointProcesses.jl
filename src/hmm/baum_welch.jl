function baum_welch_multiple_sequences(
    hmm::HMM{TransitionsType,EmissionsType}, observation_sequences; iterations
) where {TransitionsType,EmissionsType}
    S = nb_states(hmm)
    K = length(observation_sequences)
    T = [length(observation_sequences[k]) for k in 1:K]
    concat_observation_sequences = reduce(vcat, observation_sequences)

    obs_density = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    α = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    β = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    c = [Vector{Float64}(undef, T[k]) for k in 1:K]
    γ = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    ξ = [Array{Float64,3}(undef, T[k] - 1, S, S) for k in 1:K]

    logL_evolution = Float64[]
    @showprogress "Baum-Welch algorithm iterations: " for _ in 1:iterations
        logL = 0.0
        for k in 1:K
            update_obs_density!(obs_density[k], observation_sequences[k], hmm)
            logL += forward_backward!(α[k], β[k], c[k], γ[k], ξ[k], obs_density[k], hmm)
        end
        push!(logL_evolution, logL)

        concat_γ = [mapreduce(x -> x[:, s], vcat, γ) for s in 1:nb_states(hmm)]
        new_transitions = fit_mle(TransitionsType, γ, ξ)
        new_emissions = [
            fit_mle(EmissionsType, concat_observation_sequences, concat_γ[s]) for
            s in 1:nb_states(hmm)
        ]
        hmm = HMM(new_transitions, new_emissions)
    end

    return hmm, logL_evolution
end

function baum_welch_multiple_sequences_log(
    hmm::HMM{TransitionsType,EmissionsType}, observation_sequences; iterations
) where {TransitionsType,EmissionsType}
    S = nb_states(hmm)
    K = length(observation_sequences)
    T = [length(observation_sequences[k]) for k in 1:K]
    concat_observation_sequences = reduce(vcat, observation_sequences)

    obs_logdensity = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    logα = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    logβ = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    logγ = [Matrix{Float64}(undef, T[k], S) for k in 1:K]
    logξ = [Array{Float64,3}(undef, T[k] - 1, S, S) for k in 1:K]

    logL_evolution = Float64[]
    @showprogress "Baum-Welch algorithm iterations (log scale): " for _ in 1:iterations
        logL = 0.0
        for k in 1:K
            update_obs_logdensity!(obs_logdensity[k], observation_sequences[k], hmm)
            logL += forward_backward_log!(
                logα[k], logβ[k], logγ[k], logξ[k], obs_logdensity[k], hmm
            )
        end
        push!(logL_evolution, logL)

        γ = map(x -> exp.(x), logγ)
        ξ = map(x -> exp.(x), logξ)
        new_transitions = fit_mle(TransitionsType, γ, ξ)
        concat_γ = [mapreduce(x -> x[:, s], vcat, γ) for s in 1:nb_states(hmm)]
        new_emissions = [
            fit_mle(EmissionsType, concat_observation_sequences, concat_γ[s]) for
            s in 1:nb_states(hmm)
        ]
        hmm = HMM(new_transitions, new_emissions)
    end

    return hmm, logL_evolution
end

function baum_welch(hmm, observation_sequence; log=false, plot=false, kwargs...)
    if log
        hmm, logL_evolution = baum_welch_multiple_sequences_log(
            hmm, [observation_sequence]; kwargs...
        )
    else
        hmm, logL_evolution = baum_welch_multiple_sequences(
            hmm, [observation_sequence]; kwargs...
        )
    end
    if plot
        println(
            lineplot(logL_evolution; xlabel="Baum-Welch iteration", ylabel="Log-likelihood")
        )
    end
    return hmm, logL_evolution
end
