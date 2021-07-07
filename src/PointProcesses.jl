module PointProcesses

include("history.jl")
export History, push!, append!, nb_events, has_events

include("point_process.jl")
export PointProcess, ground_intensity, ground_intensity_bound, ground_intensity_bound_validity_duration, mark_distribution, mark_density, intensity

include("ogata.jl")
export simulate

include("poisson.jl")
export PoissonProcess

include("plot.jl")
export plot

end
