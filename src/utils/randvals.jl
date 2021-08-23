"""
    uniformprobvec(n)

Return a uniform probability distribution vector of size `n`.
"""
function uniformprobvec(n)
    return ones(n) ./ n
end

"""
    uniformtransmat(n)

Return a stochastic matrix of size `n` with uniform transition probability distributions.
"""
function uniformtransmat(n)
    return ones(n, n) ./ n
end

"""
    randprobvec(n)

Return a random probability distribution vector of size `n`.
"""
function randprobvec(n)
    p = rand(n)
    return p ./ sum(p)
end

"""
    randtransmat(n)

Return a stochastic matrix of size `n` with random transition probability distributions.
"""
function randtransmat(n)
    P = rand(n, n)
    return P ./ sum(P, dims = 2)
end
