"""
    AbstractPointProcess{M,L}

The common supertype of all point processes with mark type `M` and location type `L`.
"""
abstract type AbstractPointProcess{M,L} <: AbstractMeasure end

"""
    build_transform(pp)

Return a transformation object from `TransformVariables` that can turn a `Vector{Float64}` into a `NamedTuple` with fields matching those of `pp`.
"""
function build_transform end
