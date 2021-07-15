module PointProcesses

include("history.jl")
export History, push!, append!, nb_events, has_events

include("point_process.jl")
export PointProcess, ground_intensity, intensity, mark_distribution, ground_intensity_bound, ground_intensity_bound_validity_duration

include("ogata.jl")
export rand

include("poisson.jl")
export PoissonProcess, default_param

include("learning.jl")
export integrated_ground_intensity, logpdf, fit

include("plot.jl")
export plot

end
