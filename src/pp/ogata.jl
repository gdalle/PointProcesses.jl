"""
    rand(rng, pp, tmin, tmax)

Simulate a point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm[^Rasmussen_2018].
"""
function Base.rand(rng::AbstractRNG, pp::PointProcess{M}, tmin::Float64, tmax::Float64) where {M}
    h = History(Float64[], M[], tmin, tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, h, t)
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
