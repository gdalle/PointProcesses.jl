sampletype(::Dists.DiscreteUnivariateDistribution) = Int
sampletype(::Dists.DiscreteMultivariateDistribution) = Vector{Int}
sampletype(::Dists.ContinuousUnivariateDistribution) = Float64
sampletype(::Dists.ContinuousMultivariateDistribution) = Vector{Float64}
