"""
    rand(rng, pp, tmin, tmax)

Simulate a temporal point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm.
"""
function Base.rand(rng::AbstractRNG, pp::TemporalPointProcess{M}, tmin, tmax) where {M}
    h = History(times=Float64[], marks=M[], tmin=tmin, tmax=tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, h, t + eps(t))
        T = B > 0 ? rand(rng, Exponential(1 / B)) : Inf
        if T > L
            t = t + L
        elseif T <= L
            U = rand(Uniform(0, 1))
            if U < ground_intensity(pp, h, t + T) / B
                m = rand(rng, mark_distribution(pp, h, t + T))
                if t + T < tmax
                    push!(h, t + T, m)
                end
            end
            t = t + T
        end
    end
    return h
end

function Base.rand(rng::AbstractRNG, tpp::BoundedTemporalPointProcess)
    return rand(rng, tpp.pp, tpp.tmin, tpp.tmax)
end

function Base.rand(tpp::TemporalPointProcess, args...; kwargs...)
    return rand(GLOBAL_RNG, tpp, args...; kwargs...)
end
