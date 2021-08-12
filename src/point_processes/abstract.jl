"""
    AbstractPointProcess{L, M}

The common supertype of all point processes with location type `L` and mark type `M`.
"""
abstract type AbstractPointProcess{L,M} end

"""
    build_transform(pp)

Return a transformation object from `TransformVariables` that can turn a `Vector{Float64}` into a `NamedTuple` with fields matching those of `pp`.
"""
function build_transform end
