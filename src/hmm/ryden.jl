function intensity_by_state(mmpp::MarkovModulatedPoissonProcess, m=nothing)
    λ = [intensity(emission(mmpp, s), m) for s = 1:nstates(mmpp)]
    Λ = Diagonal(λ)
    return Λ
end

function renewal_density(mmpp::MarkovModulatedPoissonProcess, y, m=nothing)
    Q = rate_matrix(mmpp)
    Λ = intensity_by_state(mmpp)
    Λm = intensity_by_state(mmpp, m)
    return exp((Q-Λ)*y) * Λm
end

function renewal_void_probability(mmpp::MarkovModulatedPoissonProcess, y)
    Q = rate_matrix(mmpp)
    Λ = intensity_by_state(mmpp)
    return exp((Q-Λ)*y)
end

function stationary_distribution(mmpp::MarkovModulatedPoissonProcess)
    Q = rate_matrix(mmpp)
    Λ = intensity_by_state(mmpp)
    P = inv(Λ - Q) * Λ
    return stationary_distribution(DiscreteMarkovChain(ones(nstates(mmpp)), P))
end


function forward_backward(
    mmpp::MarkovModulatedPoissonProcess,
    h::TemporalHistory,
)
    S = size(Q, 1)
    π0 = stationary_distribution(mmpp)

    Y = diff(vcat(h.t_min, h.events))
    marks = h.marks
    n = nb_events(h)

    L = OffsetArray{Adjoint{Float64,Array{Float64,1}},1}(undef, 0:n+1)
    R = Vector{Vector{Float64}}(undef, n + 2)
    c = Vector{Float64}(undef, n + 1)
    C = Vector{Matrix{Float64}}(undef, n + 1)
    ∫ = Vector{Matrix{Float64}}(undef, n + 1)

    f_values = Matrix{Float64}[f(Y[k], Q, λ, m = marks[k]) for k = 1:n]

    L[0] = π0
    for k = 1:n+1
        if k < n + 1
            L[k] = L[k-1] * f_values[k]
        else
            L[n+1] = L[n] * F̄(h.t_max - h.events[end], Q, λ)
        end
        c[k] = sum(L[k])
        L[k] ./= c[k]
    end

    R[n+2] = ones(Float64, S)
    for k = n+1:-1:1
        if k == n + 1
            R[n+1] = F̄(h.t_max - h.events[end], Q, λ) * R[n+2]
        else
            R[k] = f_values[k] * R[k+1]
        end
        R[k] ./= c[k]
    end

    for k = 1:n+1
        C[k] = Matrix{Float64}(undef, 2S, 2S)
        C[k][1:S, 1:S] = Q - get_Λ(λ)
        C[k][S+1:2S, 1:S] .= 0.0
        C[k][S+1:2S, S+1:2S] = Q - get_Λ(λ)
        if k < n + 1
            C[k][1:S, S+1:2S] = get_Λ(λ, m = marks[k]) * R[k+1] * L[k-1]
        else
            C[k][1:S, S+1:2S] = R[k+1] * L[k-1]
        end
    end

    for k = 1:n+1
        if k < n + 1
            ∫[k] = exp(C[k] * Y[k])[1:S, S+1:2S]
        else
            ∫[n+1] = exp(C[n+1] * (h.t_max - h.events[end]))[1:S, S+1:2S]
        end
    end

    return L, R, c, ∫
end

function rate_loglikelihood(
    λ::CategoricalRates;
    α_prior::Matrix{Float64},
    β_prior::Matrix{Float64}
)
    S, M = size(λ.λ)
    LL = 0.
    for s = 1:S, m = 1:M
        LL += logpdf(Gamma(α_prior[s, m], 1/β_prior[s, m]), λ.λ[s, m])
    end
    return LL
end

function Estep(
    h::MarkedHistory{MarkType},
    Q::Matrix{Float64},
    λ::Rates,
)::NamedTuple where MarkType
    S = size(Q, 1)

    D̂ = zeros(Float64, S)
    m̂ = zeros(Float64, S, S)
    n̂ = DefaultDict{MarkType,Vector{Float64}}(zeros(Float64, S))
    LL = 0.0

    L, R, c, ∫ = forward_backward(h, Q, λ)

    events, marks = h.events, h.marks
    n = length(events)

    m̂_loc = Q .* sum(∫[k]' / c[k] for k = 1:n+1)
    D̂_loc = diag(m̂_loc) ./ diag(Q)
    m̂ += m̂_loc
    D̂ += D̂_loc

    for k = 1:n
        n̂[marks[k]] += L[k]' .* R[k+1]
    end

    LL += sum(log.(c))

    return (D̂=D̂, m̂=m̂, n̂=n̂, LL=LL)
end

function Estep!(
    histories::Vector{<:MarkedHistory},
    Q::Matrix{Float64},
    λ::Rates,
    suffstats::Vector{<:NamedTuple}
)
    for (r, h) in enumerate(histories)
        suffstats[r] = Estep(h, Q, λ)
    end
end

function Estep_incremental!(
    histories::Vector{<:MarkedHistory},
    Q::Matrix{Float64},
    λ::Rates,
    suffstats::Vector{<:NamedTuple};
    batch_size=nothing
)
    if batch_size === nothing
        batch_size = length(histories)
    end
    for r in sample(1:length(histories), batch_size, replace=false)
        suffstats[r] = Estep(histories[r], Q, λ)
    end
end

## M step

function Mstep(
    suffstats::Vector{<:NamedTuple},
    old_λ::Rates;
    α_prior::Matrix{Float64},
    β_prior::Matrix{Float64},
)
    D̂ = sum(ss[:D̂] for ss in suffstats)
    m̂ = sum(ss[:m̂] for ss in suffstats)
    n̂ = sum(ss[:n̂] for ss in suffstats)
    LL = sum(ss[:LL] for ss in suffstats) + rate_loglikelihood(old_λ, α_prior=α_prior, β_prior=β_prior)

    new_Q = m̂ ./ D̂
    new_Q[diagind(new_Q)] .= 0.0
    new_Q[diagind(new_Q)] = -sum(new_Q, dims = 2)
    S, M = size(old_λ.λ)
    new_λ = CategoricalRates(collect(
        (n̂[m][s] + α_prior[s, m] - 1) / (D̂[s] + β_prior[s, m]) for s = 1:S, m = 1:M
    ))

    return new_Q, new_λ, LL
end

function Mstep_marked(
    D̂::Vector{Float64},
    m̂::Matrix{Float64},
    n̂::DefaultDict{Vector{Int64}},
    old_λ::MarkedRates,
)::Tuple{Matrix{Float64},MarkedRates}
    new_Q = m̂ ./ D̂
    new_Q[diagind(new_Q)] .= 0.0
    new_Q[diagind(new_Q)] = -sum(new_Q, dims = 2)

    S, D = size(old_λ.p)
    mark_dims = collect(length(old_λ.p[1, d]) for d = 1:D)
    new_p = Matrix{Vector{Float64}}(undef, S, D)
    for s = 1:S, d = 1:D
        new_p[s, d] = Vector{Float64}(undef, mark_dims[d])
    end

    all_n̂ = sum(values(n̂))
    for d = 1:D
        for val = 1:mark_dims[d]
            right_n̂ =
                sum([n̂[m] for m in filter(x -> (x[d] == val), keys(n̂))])
            for s = 1:S
                new_p[s, d][val] = right_n̂[s] / all_n̂[s]
            end
        end
    end
    new_λ = MarkedRates(sum(values(n̂)) ./ D̂, new_p)
    return new_Q, new_λ
end
#
# function Mstep_opt(
#     D̂::Vector{Float64},
#     m̂::Matrix{Float64},
#     n̂::Matrix{Float64},
# )::Tuple{Matrix{Float64},Matrix{Float64}}
#     new_Q = m̂ ./ D̂
#     new_Q[diagind(new_Q)] .= 0.0
#     new_Q[diagind(new_Q)] = -sum(new_Q, dims = 2)
#
#     new_λ = Variable(S, M)
#     problem = maximize(
#         dot(n̂ .+ 1e-5, log(new_λ)) - dot(D̂ .+ 1e-5, sum(new_λ, dims = 2)),
#         [new_λ > 0],
#     )
#     solve!(problem, SCS.Optimizer(verbose = false))
#     return new_Q, new_λ.value
# end
#
# function Mstep_split(
#     D̂::Vector{Float64},
#     m̂::Matrix{Float64},
#     n̂::Matrix{Float64},
# )::Tuple{Matrix{Float64},Matrix{Float64}}
#     new_Q = m̂ ./ D̂
#     new_Q[diagind(new_Q)] .= 0.0
#     new_Q[diagind(new_Q)] = -sum(new_Q, dims = 2)
#
#     new_λ = Variable(S, M)
#     λ1 = Variable(2, M)
#     λ2 = Variable(S ÷ 2, M)
#     objective =
#         dot(n̂ .+ 1e-5, log(new_λ)) - dot(D̂ .+ 1e-5, sum(new_λ, dims = 2))
#     constraints = Constraint[λ1>0, λ2>0]
#     for s = 0:S-1
#         push!(constraints, new_λ[s+1, :] == λ1[(s%2)+1, :] + λ2[(s÷2)+1, :])
#     end
#     problem = maximize(objective, constraints)
#     solve!(problem, SCS.Optimizer(verbose = false))
#     return new_Q, new_λ.value
# end

## EM algorithm

function em_iteration!(
    histories::Vector{<:MarkedHistory},
    Q::Matrix{Float64},
    λ::Rates,
    suffstats::Vector{<:NamedTuple};
    α_prior::Matrix{Float64},
    β_prior::Matrix{Float64},
    batch_size=nothing
)::Tuple{Matrix{Float64},Rates,Float64}

    Estep_incremental!(histories, Q, λ, suffstats, batch_size = batch_size)
    new_Q, new_λ, LL = Mstep(suffstats, λ, α_prior = α_prior, β_prior = β_prior)

    return new_Q, new_λ, LL
end

function em(
    histories::Vector{<:MarkedHistory},
    Q::Matrix{Float64},
    λ::Rates;
    α_prior::Matrix{Float64},
    β_prior::Matrix{Float64},
    iterations::Int64 = 100,
    threshold = 1e-2,
    batch_size = nothing,
)::Tuple{Matrix{Float64},Rates,Vector{Float64}}

    suffstats = Vector{NamedTuple}(undef, length(histories))
    for (r, h) in enumerate(histories)
        D̂ = zeros(Float64, S)
        m̂ = zeros(Float64, S, S)
        n̂ = DefaultDict{Any,Vector{Float64}}(zeros(Float64, S))
        LL = 0.0
        suffstats[r] = (D̂=D̂, m̂=m̂, n̂=n̂, LL=LL)
    end

    prog = ProgressThresh(threshold, "Minimizing loglikelihood change:")
    loglikelihood_evolution = Float64[]

    for it = 1:iterations
        try
            local_batch_size = (it >= 2) ? batch_size : nothing
            new_Q, new_λ, LL = em_iteration!(
                histories,
                Q,
                λ,
                suffstats,
                α_prior = α_prior,
                β_prior = β_prior,
                batch_size = local_batch_size
            )
            conv_test = it >= 2 ? LL - loglikelihood_evolution[end] : Inf
            push!(loglikelihood_evolution, LL)
            update!(prog, conv_test)
            Q, λ = new_Q, new_λ
            if conv_test < threshold
                break
            end
        catch e
            if isa(e, InterruptException)
                break
            else
                rethrow(e)
            end
        end

    end

    # @assert minimum(diff(loglikelihood_evolution)) > 0.0
    plot_loglikelihood(loglikelihood_evolution)

    return Q, λ, loglikelihood_evolution
end

function learn_MMPP(
    histories::Vector{<:MarkedHistory},
    S::Int64;
    Q_init,
    λ_init,
    α_prior::Matrix{Float64},
    β_prior::Matrix{Float64},
    iterations::Int64 = 100,
    threshold = 1e-2,
    batch_size = nothing,
)::Tuple{Matrix{Float64},Rates,Vector{Float64}}
    Q_hat, λ_hat, loglikelihood_evolution = em(
        histories,
        Q_init,
        λ_init,
        α_prior = α_prior,
        β_prior = β_prior,
        iterations = iterations,
        threshold = threshold,
        batch_size = batch_size,
    )
    return Q_hat, λ_hat, loglikelihood_evolution
end

## State inference

function infer_state(
    h::MarkedHistory{Int64},
    Q::Matrix{Float64},
    λ::Rates;
    times::Vector{Float64},
    plot::Bool = true,
)
    events, marks = h.events, h.marks
    n = nb_events(h)
    completed_events = OffsetArray(vcat(h.t_min, events, h.t_max), 0:n+1)
    S = size(Q, 1)
    L, R, c, ∫ = forward_backward(h, Q, λ)
    ŝ = Matrix{Float64}(undef, length(times), S)
    for (nt, t) in enumerate(times)
        k = searchsortedfirst(completed_events, t)
        t1, t2 = completed_events[k-1], completed_events[k]
        @assert t1 < t < t2
        transition_before = L[k-1] * F̄(t - t1, Q, λ)
        transition_after = k < n + 1 ? f(t2 - t, Q, λ, m = marks[k]) * R[k+1] :
            F̄(t2 - t, Q, λ) * R[k+1]
        ŝ[nt, :] = transition_before' .* transition_after
        ŝ[nt, :] ./= sum(ŝ[nt, :])
    end
    ŝ_max = [argmax(ŝ[nt, :]) for nt in eachindex(times)]
    if plot
        plot_state_inference(times, ŝ, ŝ_max)
    end
    return ŝ
end

function infer_last_state(h::MarkedHistory{Int64}, Q::Matrix{Float64}, λ::Rates)
    ŝ = infer_state(h, Q, λ, times = [h.t_max - ε], plot = false)
    return ŝ[1, :]
end

function predict_next_event(
    h::MarkedHistory{Int64},
    Q::Matrix{Float64},
    λ::Rates,
    m::Int64,
)
    ŝ = infer_last_state(h, Q, λ)
    return ŝ' * inv(get_Λ(λ, m = m) - Q)^2 * get_Λ(λ, m = m) * ones(size(Q, 1))
end
