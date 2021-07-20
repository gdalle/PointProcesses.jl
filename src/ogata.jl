function Distributions.rand(pp::PointProcess{M}, tmin, tmax) where {M}
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
