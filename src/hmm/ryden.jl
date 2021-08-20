function intensity_by_state(mmpp::MMPP, m = nothing)
    λ = [intensity(emission(mmpp, s), m) for s = 1:nstates(mmpp)]
    Λ = Diagonal(λ)
    return Λ
end

function renewal_density(mmpp::MMPP, y, m = nothing)
    Q = rate_matrix(mmpp)
    Λ = intensity_by_state(mmpp)
    Λm = intensity_by_state(mmpp, m)
    return exp((Q - Λ) * y) * Λm
end

function renewal_void_probability(mmpp::MMPP, y)
    Q = rate_matrix(mmpp)
    Λ = intensity_by_state(mmpp)
    return exp((Q - Λ) * y)
end

function stationary_distribution(mmpp::MMPP)
    Q = rate_matrix(mmpp)
    Λ = intensity_by_state(mmpp)
    P = inv(Λ - Q) * Λ
    return stationary_distribution(DiscreteMarkovChain(ones(nstates(mmpp)), P))
end

function forward_backward(mmpp::MMPP, h::TemporalHistory)
    S = nstates(mmpp)
    π0 = stationary_distribution(mmpp)

    y = diff(vcat(min_time(h), event_times(h), max_time(h)))
    marks = event_marks(h)
    n = nb_events(h)

    Λ = intensity_by_state(mmpp)
    f = [renewal_density(mmpp, y[k], marks[k]) for k = 1:n]
    F̄ = renewal_void_probability(mmpp, y[n+1])

    L = OffsetVector{Float64}(undef, 0:n+1)
    c = Vector{Float64}(undef, n + 1)
    L[0] = π0
    for k = 1:n+1
        L[k] = k < n + 1 ? L[k-1] * f[k] : L[n] * F̄
        c[k] = sum(L[k])
        L[k] ./= c[k]
    end

    R = Vector{Vector{Float64}}(undef, n + 2)
    R[n+2] = ones(S)
    for k = n+1:-1:1
        R[n+1] = k == n + 1 ? F̄ * R[n+2] : f[k] * R[k+1]
        R[k] ./= c[k]
    end

    C = Matrix{Float64}(undef, 2S, 2S)
    C[1:S, 1:S] = Q - Λ
    C[S+1:2S, 1:S] .= 0.0
    C[S+1:2S, S+1:2S] = Q - Λ

    ∫ = Vector{Matrix{Float64}}(undef, n + 1)
    for k = 1:n+1
        if k < n + 1
            Λm = intensity_by_state(mmpp, marks[k])
            C[1:S, S+1:2S] = Λm * R[k+1] * L[k-1]
        else
            C[1:S, S+1:2S] = R[k+1] * L[k-1]
        end
        ∫[k] = exp(C[k] * y[k])[1:S, S+1:2S]
    end

    Q = rate_matrix(mmpp)

    m̂ = Q .* sum(∫[k]' / c[k] for k = 1:n+1)
    D̂ = diag(m̂) ./ diag(Q)
    n̂ = [DefaultDict{M,Float64}(0.0) for s = 1:S]
    for k = 1:n
        local_n̂ = L[k]' .* R[k+1]
        for s = 1:S
            n̂[s][marks[k]] += local_n̂[s]
        end
    end

    logL = sum(log.(c))

    return m̂, n̂, D̂, logL
end

function ryden(mmpp::MMPP{M,Tr,Em}, h::TemporalHistory{M}; iterations) where {M,Tr,Em}
    logL_evolution = Float64[]
    for _ = 1:iterations
        m̂, n̂, D̂, logL = forward_backward(mmpp, h)
        push!(logL_evolution, logL)
        new_transitions = fit(Tr, m̂ = m̂, D̂ = D̂)
        new_emissions = [fit(Em, n̂ = n̂[s], D̂ = D̂[s]) for s = 1:nstates(hmm)]  # TODO: @view
        mmpp = MMPP(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end
