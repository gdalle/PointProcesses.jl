function uniformprobvec(n)
    return ones(n) ./ n
end

function uniformtransmat(n)
    return ones(n, n) ./ n
end

function randprobvec(n)
    p = rand(n)
    return p ./ sum(p)
end

function randtransmat(n)
    P = rand(n, n)
    return P ./ sum(P, dims=2)
end
