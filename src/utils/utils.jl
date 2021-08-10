"""
    logsumexp(x)

Compute the `logsumexp` function on any array of real numbers without numerical instability.
"""
function logsumexp(x::AbstractArray{R}) where {R<:Real}
    n = length(x)
    if n == 0
        return typemin(R)
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
        if y > typemin(y)
            return false
        end
    end
    return true
end

function all_plus_inf(x)
    for y in x
        if y < typemax(y)
            return false
        end
    end
    return true
end

function all_zeros(x)
    for y in x
        if 1 / abs(x) < typemax(x)
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
