Base.eltype(::Type{<:TemporalPointProcess{M}}) where {M} = TemporalHistory{M}

"""
    rand(rng, pp, tmin, tmax)

Simulate a temporal point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm[^Ogata_1981].

[^Ogata_1981]: Ogata, Y. (1981), “On Lewis’ simulation method for point processes,” IEEE Transactions on Information Theory, 27, 23–31. https://doi.org/10.1109/TIT.1981.1056305.
"""
function Base.rand(rng::AbstractRNG, pp::TemporalPointProcess{M}, tmin, tmax) where {M}
    h = TemporalHistory(times = Float64[], marks = M[], tmin = tmin, tmax = tmax)
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

Base.rand(rng::AbstractRNG, tpp::BoundedTemporalPointProcess) =
    rand(rng, tpp.pp, tpp.tmin, tpp.tmax)

Base.rand(tpp::TemporalPointProcess, args...) = rand(Random.GLOBAL_RNG, tpp, args...)
