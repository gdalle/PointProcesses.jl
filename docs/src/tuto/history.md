# Event history

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

To analyze point processes, we need a way to store their realizations. This is the purpose of the [`AbstractHistory`](@ref) subtypes, of which only [`TemporalHistory`](@ref) is implemented for now.

```jldoctest
julia> history = TemporalHistory([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0);

julia> duration(history)
2.0

julia> nb_events(history)
3

julia> nb_events(history, 1.0, 2.0)
1

julia> has_events(history)
true

julia> has_events(history, 1.5, 2.0)
false

julia> push!(history, 1.7, "d")

julia> has_events(history, 1.5, 2.0)
true
```