# Event history

```@setup
ENV["GKSwstype"] = "100"
```

To analyze point processes, we need a way to store their realizations. This is the purpose of the [`AbstractHistory`](@ref) subtypes, of which only [`TemporalHistory`](@ref) is implemented for now.

```@repl
using PointProcesses

history = TemporalHistory([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0);
duration(history)
nb_events(history)
nb_events(history, 1.0, 2.0)
has_events(history)
has_events(history, 1.5, 2.0)
push!(history, 1.7, "d")
has_events(history, 1.5, 2.0)
```