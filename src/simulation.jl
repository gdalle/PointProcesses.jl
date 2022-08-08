"""
    rand(rng, pp, tmin, tmax)

Simulate a temporal point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm.
"""
function Base.rand(
    rng::AbstractRNG, pp::AbstractPointProcess{M}, tmin::T, tmax::T
) where {M,T<:Real}
    h = History(; times=T[], marks=M[], tmin=tmin, tmax=tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, h, t + eps(t))
        τ = B > 0 ? rand(rng, Exponential(inv(B))) : Inf
        if τ > L
            t = t + L
        elseif τ <= L
            U = rand(rng, T)
            if U < ground_intensity(pp, h, t + τ) / B
                m = rand(rng, mark_distribution(pp, h, t + τ))
                if t + τ < tmax
                    push!(h, t + τ, m)
                end
            end
            t = t + τ
        end
    end
    return h
end

function Base.rand(pp::AbstractPointProcess, args...; kwargs...)
    return rand(GLOBAL_RNG, pp, args...; kwargs...)
end
