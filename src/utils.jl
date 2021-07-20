"""
    logsumexp(x)

Compute the `logsumexp` function on any array of real numbers without numerical instability.
"""
function logsumexp(x::AbstractArray{<:Real})
    n = length(x)
    if n == 0
        return -Inf
    else
        m = maximum(x)
        s = exp(x[1] - m)
        for i = 2:n
            s += exp(x[i] - m)
        end
        return m + log(s)
    end
end

# Overflow checks

function all_minus_inf(x)
    for y in x
        if y > -Inf
            return false
        end
    end
    return true
end

function all_plus_inf(x)
    for y in x
        if y < +Inf
            return false
        end
    end
    return true
end

function all_zeros(x)
    for y in x
        if 1/abs(x) < Inf
            return false
        end
    end
    return true
end

function all_nan(x)
    for y in x
        if !isnan(y)
            return false
        end
    end
    return true
end

function Plots.scatter(history::History{M}) where {M}
    scatter(
        history.times,
        history.marks,
        xlim = (history.tmin, history.tmax),
        title = "Event history",
        xlabel = "Time",
        ylabel = "Marks",
        label = nothing,
    )
end
