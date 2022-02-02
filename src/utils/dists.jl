sampletype(::DiscreteUnivariateDistribution) = Int
sampletype(::DiscreteMultivariateDistribution) = Vector{Int}
sampletype(::ContinuousUnivariateDistribution) = Float64
sampletype(::ContinuousMultivariateDistribution) = Vector{Float64}
