using Distributions
import Distributions: rand

function extend_ogata!(
    pointprocess::PointProcess{MarkType},
    history::History{MarkType},
    tmax,
) where {MarkType}
    t = history.tmax
    while t < tmax
        M = ground_intensity_bound(pointprocess, history, t)
        L = ground_intensity_bound_validity_duration(pointprocess, history, t)
        T = M > 0 ? rand(Exponential(1 / M)) : Inf
        if T > L
            t = t + L
        elseif T <= L
            U = rand(Uniform(0, 1))
            if U < ground_intensity(pointprocess, history, t + T) / M
                m = rand(mark_distribution(pointprocess, history, t + T))
                if t + T < tmax
                    push!(history, t + T, m)
                end
            end
            t = t + T
        end
    end
    history.tmax = tmax
    return nothing
end

function simulate(pointprocess::PointProcess{MarkType}; tmin, tmax) where {MarkType}
    history = History(Float64[], MarkType[], tmin, tmin)
    extend_ogata!(pointprocess, history, tmax)
    return history
end
