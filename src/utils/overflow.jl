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
        if 1 / abs(y) < typemax(y)
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
