function baum_welch!(
    α,
    β,
    c,
    γ,
    ξ,
    obs_density::AbstractMatrix,
    hmm::HMM{TransitionsType,EmissionsType},
    observations::AbstractVector;
    iterations,
) where {TransitionsType,EmissionsType}
    logL_evolution = Float64[]
    @showprogress "Baum-Welch algorithm iterations: " for _ in 1:iterations
        update_obs_density!(obs_density, hmm, observations)
        logL = forward_backward!(α, β, c, γ, ξ, obs_density, hmm)
        push!(logL_evolution, logL)
        new_transitions = fit_mle(TransitionsType; γ=γ, ξ=ξ)
        new_emissions = [
            fit_mle(EmissionsType, observations, γ[:, s]) for s in 1:nb_states(hmm)
        ]  # TODO: @view
        hmm = HMM(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end

function baum_welch_log!(
    logα,
    logβ,
    logγ,
    logξ,
    obs_logdensity::AbstractMatrix,
    hmm::HMM{TransitionsType,EmissionsType},
    observations::AbstractVector;
    iterations,
) where {TransitionsType,EmissionsType}
    logL_evolution = Float64[]
    @showprogress "Baum-Welch algorithm iterations (log scale): " for _ in 1:iterations
        update_obs_logdensity!(obs_logdensity, hmm, observations)
        logL = forward_backward_log!(logα, logβ, logγ, logξ, obs_logdensity, hmm)
        push!(logL_evolution, logL)
        new_transitions = fit_mle(TransitionsType; γ=exp.(logγ), ξ=exp.(logξ))
        new_emissions = [
            fit_mle(EmissionsType, observations, exp.(logγ[:, s])) for s in 1:nb_states(hmm)
        ]  # TODO: @view
        hmm = HMM(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end

## Non mutating version

function baum_welch(
    hmm::HMM, observations::AbstractVector; iterations, log=false, plot=false
)
    T, S = length(observations), nb_states(hmm)
    if log
        logα = Matrix{Float64}(undef, T, S)
        logβ = Matrix{Float64}(undef, T, S)
        logγ = Matrix{Float64}(undef, T, S)
        logξ = Array{Float64,3}(undef, T - 1, S, S)
        obs_logdensity = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution = baum_welch_log!(
            logα, logβ, logγ, logξ, obs_logdensity, hmm, observations; iterations=iterations
        )
    else
        α = Matrix{Float64}(undef, T, S)
        β = Matrix{Float64}(undef, T, S)
        c = Vector{Float64}(undef, T)
        γ = Matrix{Float64}(undef, T, S)
        ξ = Array{Float64,3}(undef, T - 1, S, S)
        obs_density = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution = baum_welch!(
            α, β, c, γ, ξ, obs_density, hmm, observations; iterations=iterations
        )
    end
    if plot
        println(
            lineplot(logL_evolution; xlabel="Baum-Welch iteration", ylabel="Loglikelihood")
        )
    end
    return hmm_est, logL_evolution
end
