function Base.rand(rng::AbstractRNG, pp::PointProcess{M}, tmin::Real, tmax::Real) where {M}
    h = History(Float64[], M[], tmin, tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, h, t)
        T = B > 0 ? rand(Exponential(1 / B)) : Inf
        if T > L
            t = t + L
        elseif T <= L
            U = rand(Uniform(0, 1))
            if U < ground_intensity(pp, h, t + T) / B
                m = rand(mark_distribution(pp, h, t + T))
                if t + T < tmax
                    push!(h, t + T, m)
                end
            end
            t = t + T
        end
    end
    return h
end

"""
    rand(pp, tmin, tmax)

Simulate a point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm[^Rasmussen_2018].

# Examples

```jldoctest
julia> pp = MultivariatePoissonProcess([0., 1., 2.]);

julia> h = rand(pp, 0., 1000.);
```
"""
Base.rand(pp::PointProcess, tmin::Real, tmax::Real) = rand(GLOBAL_RNG, pp, tmin, tmax)
