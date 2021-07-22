Base.eltype(::Type{<:PointProcess{M}}) where {M} = History{M}

"""
    rand(rng, pp, tmin, tmax)

Simulate a point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm[^Ogata_1981].

[^Ogata_1981]: Ogata, Y. (1981), “On Lewis’ simulation method for point processes,” IEEE Transactions on Information Theory, 27, 23–31. https://doi.org/10.1109/TIT.1981.1056305.
"""
function Base.rand(rng::AbstractRNG, pp::PointProcess{M}, tmin::Float64, tmax::Float64) where {M}
    h = History(Float64[], M[], tmin, tmax)
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

Base.rand(rng::AbstractRNG, tpp::TimedPointProcess) = rand(rng, tpp.pp, tpp.tmin, tpp.tmax)

Base.rand(pp::PointProcess, args...) = rand(GLOBAL_RNG, pp, args...)
