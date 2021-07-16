using Distributions
import Distributions: rand

function rand(pp::PointProcess{M}, tmin, tmax) where {M}
    history = History(Float64[], M[], tmin, tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, history, t)
        T = B > 0 ? rand(Exponential(1 / B)) : Inf
        if T > L
            t = t + L
        elseif T <= L
            U = rand(Uniform(0, 1))
            if U < ground_intensity(pp, history, t + T) / B
                m = rand(mark_distribution(pp, history, t + T))
                if t + T < tmax
                    push!(history, t + T, m)
                end
            end
            t = t + T
        end
    end
    return history
end