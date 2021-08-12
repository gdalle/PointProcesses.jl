var documenterSearchIndex = {"docs":
[{"location":"history/#Event-history","page":"Event history","title":"Event history","text":"","category":"section"},{"location":"history/","page":"Event history","title":"Event history","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"history/abstract.jl\"]","category":"page"},{"location":"history/#PointProcesses.AbstractHistory","page":"Event history","title":"PointProcesses.AbstractHistory","text":"AbstractHistory{L, M}\n\nAbstract supertype for event histories with locations of type L and marks of type M.\n\n\n\n\n\n","category":"type"},{"location":"history/","page":"Event history","title":"Event history","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"history/temporal.jl\"]","category":"page"},{"location":"history/#PointProcesses.TemporalHistory","page":"Event history","title":"PointProcesses.TemporalHistory","text":"TemporalHistory{M}\n\nLinear event histories with marks of type M.\n\nFields\n\ntimes::Vector{Float64}: vector of event times\nmarks::Vector{M}: vector of event marks\ntmin::Float64: start time\ntmax::Float64: end time\n\n\n\n\n\n","category":"type"},{"location":"history/#Base.push!-Tuple{TemporalHistory, Any, Any}","page":"Event history","title":"Base.push!","text":"push!(h::TemporalHistory{M}, t::Float64, m::M)\n\nAdd event (t, m) at the end of history h.\n\n\n\n\n\n","category":"method"},{"location":"history/#PointProcesses.duration-Tuple{TemporalHistory}","page":"Event history","title":"PointProcesses.duration","text":"duration(h::TemporalHistory)\n\nCompute the difference h.tmax - h.tmin.\n\n\n\n\n\n","category":"method"},{"location":"history/#PointProcesses.has_events","page":"Event history","title":"PointProcesses.has_events","text":"has_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)\n\nCheck the presence of events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"function"},{"location":"history/#PointProcesses.nb_events","page":"Event history","title":"PointProcesses.nb_events","text":"nb_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)\n\nCount events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"function"},{"location":"history/#PointProcesses.time_change-Tuple{TemporalHistory, Any}","page":"Event history","title":"PointProcesses.time_change","text":"time_change(h, Λ)\n\nApply the time rescaling t mapsto Lambda(t) to history h.\n\n\n\n\n\n","category":"method"},{"location":"models/#Point-process-models","page":"Built-in models","title":"Point process models","text":"","category":"section"},{"location":"models/","page":"Built-in models","title":"Built-in models","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"models/poisson.jl\"]","category":"page"},{"location":"models/","page":"Built-in models","title":"Built-in models","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"models/temporal_poisson.jl\"]","category":"page"},{"location":"models/#PointProcesses.TemporalPoissonProcess","page":"Built-in models","title":"PointProcesses.TemporalPoissonProcess","text":"TemporalPoissonProcess{R}\n\nHomogeneous temporal multivariate Poisson process.\n\nFields\n\nλ::Vector{R}: event rates.\n\n\n\n\n\n","category":"type"},{"location":"models/","page":"Built-in models","title":"Built-in models","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"models/temporal_hawkes.jl\"]","category":"page"},{"location":"models/#PointProcesses.TemporalHawkesProcess","page":"Built-in models","title":"PointProcesses.TemporalHawkesProcess","text":"TemporalHawkesProcess{R}\n\nMultivariate temporal Hawkes process.\n\nFields\n\nλ::Vector{R}: base event rates\nα::Matrix{R}: excitation amplitudes\nβ::Matrix{R}: excitation decays\n\n\n\n\n\n","category":"type"},{"location":"models/","page":"Built-in models","title":"Built-in models","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"models/cox.jl\"]","category":"page"},{"location":"point_processes/#Point-processes","page":"General functions","title":"Point processes","text":"","category":"section"},{"location":"point_processes/","page":"General functions","title":"General functions","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"point_processes/abstract.jl\"]","category":"page"},{"location":"point_processes/#PointProcesses.AbstractPointProcess","page":"General functions","title":"PointProcesses.AbstractPointProcess","text":"AbstractPointProcess{L, M}\n\nThe common supertype of all point processes with location type L and mark type M.\n\n\n\n\n\n","category":"type"},{"location":"point_processes/#Temporal-point-processes","page":"General functions","title":"Temporal point processes","text":"","category":"section"},{"location":"point_processes/","page":"General functions","title":"General functions","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"point_processes/temporal.jl\", \"point_processes/temporal_multivariate.jl\"]","category":"page"},{"location":"point_processes/#PointProcesses.BoundedTemporalPointProcess","page":"General functions","title":"PointProcesses.BoundedTemporalPointProcess","text":"BoundedTemporalPointProcess{L<:Real, M}\n\nStore a temporal point process P with pre-defined start and end times.\n\n\n\n\n\n","category":"type"},{"location":"point_processes/#PointProcesses.TemporalPointProcess","page":"General functions","title":"PointProcesses.TemporalPointProcess","text":"TemporalPointProcess{L<:Real, M}\n\nThe common supertype of all temporal point processes (i.e. point processes on the real line) with mark type M.\n\n\n\n\n\n","category":"type"},{"location":"point_processes/#PointProcesses.MultivariateTemporalPointProcess","page":"General functions","title":"PointProcesses.MultivariateTemporalPointProcess","text":"MultivariateTemporalPointProcess\n\nThe common supertype of all temporal point processes with integer marks.\n\n\n\n\n\n","category":"type"},{"location":"point_processes/#PointProcesses.all_marks-Tuple{MultivariateTemporalPointProcess}","page":"General functions","title":"PointProcesses.all_marks","text":"all_marks(pp)\n\nEnumerate possible marks for a multivariate temporal point process.\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#Intensity-functions","page":"General functions","title":"Intensity functions","text":"","category":"section"},{"location":"point_processes/","page":"General functions","title":"General functions","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"point_processes/intensity.jl\"]","category":"page"},{"location":"point_processes/#PointProcesses.ground_intensity-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}, Any}} where M","page":"General functions","title":"PointProcesses.ground_intensity","text":"ground_intensity(pp, h, t)\n\nCompute the ground intensity for a temporal point process pp applied to history h at time t.\n\nThe ground intensity quantifies the instantaneous risk of an event with any mark occurring at time t[Rasmussen_2018]. It can be expressed as\n\n    lambda_g(th) = sum_m in mathcalM lambda(t mh)\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#PointProcesses.ground_intensity_bound-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}, Any}} where M","page":"General functions","title":"PointProcesses.ground_intensity_bound","text":"ground_intensity_bound(pp, θ, h, t)\n\nCompute a local upper bound on the ground intensity for a temporal point process pp applied to history h at time t[Rasmussen_2018].\n\nReturn a tuple of the form (B L) satisfying\n\n    forall u in t t+L) quad lambda_g(th) leq B\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#PointProcesses.intensity-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}, Any, M}} where M","page":"General functions","title":"PointProcesses.intensity","text":"intensity(pp, h, t, m)\n\nCompute the conditional intensity for a temporal point process pp applied to history h and event (t, m).\n\nThe conditional intensity function lambda(t m  h) quantifies the instantaneous risk  of an event with mark m occurring at time t[Rasmussen_2018].\n\n[Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#PointProcesses.mark_distribution-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}, Any}} where M","page":"General functions","title":"PointProcesses.mark_distribution","text":"mark_distribution(pp, h, t)\n\nCompute the distribution of marks for a temporal point process pp knowing that an event takes place at time t after history h.\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#Simulation","page":"General functions","title":"Simulation","text":"","category":"section"},{"location":"point_processes/","page":"General functions","title":"General functions","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"point_processes/ogata.jl\"]","category":"page"},{"location":"point_processes/#Base.rand-Union{Tuple{M}, Tuple{Random.AbstractRNG, TemporalPointProcess{M}, Any, Any}} where M","page":"General functions","title":"Base.rand","text":"rand(rng, pp, tmin, tmax)\n\nSimulate a temporal point process pp on interval [tmin, tmax) using Ogata's algorithm[Ogata_1981].\n\n[Ogata_1981]: Ogata, Y. (1981), “On Lewis’ simulation method for point processes,” IEEE Transactions on Information Theory, 27, 23–31. https://doi.org/10.1109/TIT.1981.1056305.\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#Learning","page":"General functions","title":"Learning","text":"","category":"section"},{"location":"point_processes/","page":"General functions","title":"General functions","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"point_processes/learning.jl\"]","category":"page"},{"location":"point_processes/#Distributions.logpdf-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}}} where M","page":"General functions","title":"Distributions.logpdf","text":"logpdf(pp, h)\n\nCompute the log probability density function for a temporal point process pp applied to history h:\n\n    log f(h) = sum_i log lambda(t_i  h_i) - Lambda(h)\n\nThe default method uses a loop over events combined with integrated_ground_intensity, but it should be reimplemented for specific processes if faster computation is possible.\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#PointProcesses.check_residuals-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}}} where M","page":"General functions","title":"PointProcesses.check_residuals","text":"check_residuals(pp, h)\n\nCheck whether the point process pp is a good fit for history h by applying Ogata's time rescaling method: if (t_i)_i is a temporal point process with intensity lambda(t), then (Lambda(t_i))_i is a standard temporal Poisson process.\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#PointProcesses.integrated_ground_intensity-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}, Any}} where M","page":"General functions","title":"PointProcesses.integrated_ground_intensity","text":"integrated_ground_intensity(pp, h[, t])\n\nCompute the integrated ground intensity (or compensator) Lambda(th) for a temporal point process pp applied to history h:\n\n    Lambda(h) = int lambda_g(th) mathrmdt\n\nThe default method uses Quadrature.jl for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.\n\n\n\n\n\n","category":"method"},{"location":"point_processes/#StatsBase.fit-Union{Tuple{PP}, Tuple{M}, Tuple{PP, TemporalHistory{M}}} where {M, PP<:TemporalPointProcess{M}}","page":"General functions","title":"StatsBase.fit","text":"fit(pp0, h)\n\nCompute the optimal parameter for a temporal point process of type typeof(pp0) on history h using maximum likelihood:\n\n    hattheta = mathrmargmax f_theta(h) theta in Theta\n\nThe default method uses GalacticOptim.jl for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.\n\n\n\n\n\n","category":"method"},{"location":"hmm/#Hidden-Markov-models","page":"Hidden Markov models","title":"Hidden Markov models","text":"","category":"section"},{"location":"hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"hmm/hmm.jl\", \"hmm/forward_backward.jl\", \"hmm/baum_welch.jl\"]","category":"page"},{"location":"hmm/#PointProcesses.HiddenMarkovModel","page":"Hidden Markov models","title":"PointProcesses.HiddenMarkovModel","text":"HiddenMarkovModel{Tr, Em}\n\nHidden Markov Model with arbitrary transition model (of type Tr) and emission distributions (of type Em).\n\nFields\n\ntransitions::Tr: state evolution process.\nemissions::Vector{Em}: one emission distribution per state.\n\n\n\n\n\n","category":"type"},{"location":"markov/#Markov-chains","page":"Markov chains","title":"Markov chains","text":"","category":"section"},{"location":"markov/","page":"Markov chains","title":"Markov chains","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"markov/abstract.jl\"]","category":"page"},{"location":"markov/","page":"Markov chains","title":"Markov chains","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"markov/discrete_time.jl\"]","category":"page"},{"location":"markov/#PointProcesses.DiscreteMarkovChain","page":"Markov chains","title":"PointProcesses.DiscreteMarkovChain","text":"DiscreteMarkovChain{R<:Real}\n\nDiscrete-time Markov chain.\n\nFields\n\nπ0::Vector{R}: initial state distribution\nP::Matrix{R}: state transition matrix.\n\n\n\n\n\n","category":"type"},{"location":"markov/#PointProcesses.DiscreteMarkovChainStats","page":"Markov chains","title":"PointProcesses.DiscreteMarkovChainStats","text":"ContinuousMarkovChainStats\n\nSufficient statistics for the likelihood of a ContinuousMarkovChain.\n\n\n\n\n\n","category":"type"},{"location":"markov/","page":"Markov chains","title":"Markov chains","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"markov/continuous_time.jl\"]","category":"page"},{"location":"markov/#PointProcesses.ContinuousMarkovChain","page":"Markov chains","title":"PointProcesses.ContinuousMarkovChain","text":"ContinuousMarkovChain{R<:Real}\n\nContinuous-time Markov chain.\n\nFields\n\nπ0::Vector{R}: initial state distribution\nQ::Matrix{R}: rate matrix.\n\n\n\n\n\n","category":"type"},{"location":"markov/#PointProcesses.ContinuousMarkovChainStats","page":"Markov chains","title":"PointProcesses.ContinuousMarkovChainStats","text":"ContinuousMarkovChainStats\n\nSufficient statistics for the likelihood of a ContinuousMarkovChain.\n\n\n\n\n\n","category":"type"},{"location":"#PointProcesses","page":"Home","title":"PointProcesses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation of PointProcesses.","category":"page"},{"location":"","page":"Home","title":"Home","text":"PointProcesses","category":"page"},{"location":"#PointProcesses","page":"Home","title":"PointProcesses","text":"A package for point process modeling, simulation and inference.\n\n\n\n\n\n","category":"module"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To get started, take a look at the Tutorial.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"DocTestSetup = quote\n    using PointProcesses\nend","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In this tutorial we demonstrate the main features of PointProcesses.jl.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> using Random\n\njulia> Random.seed!(63);","category":"page"},{"location":"tutorial/#Working-with-event-histories","page":"Tutorial","title":"Working with event histories","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To analyze point processes, we need a way to store their realizations. This is the purpose of the AbstractHistory subtypes, of which only TemporalHistory is implemented for now.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> history = TemporalHistory([0.2, 0.8, 1.1], [\"a\", \"b\", \"c\"], 0.0, 2.0);\n\njulia> duration(history)\n2.0\n\njulia> nb_events(history)\n3\n\njulia> nb_events(history, 1.0, 2.0)\n1\n\njulia> has_events(history)\ntrue\n\njulia> has_events(history, 1.5, 2.0)\nfalse\n\njulia> push!(history, 1.7, \"d\")\n\njulia> has_events(history, 1.5, 2.0)\ntrue","category":"page"},{"location":"tutorial/#Working-with-Markov-processes","page":"Tutorial","title":"Working with Markov processes","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Some point processes are based on underlying Markov processes, which is why we provide a basic implementation for them.","category":"page"},{"location":"tutorial/#Discrete-time-Markov-chains","page":"Tutorial","title":"Discrete time Markov chains","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> dmc = DiscreteMarkovChain(π0 = [0.3, 0.7], P = [0.9 0.1; 0.2 0.8]);\n\njulia> states = rand(dmc, 1000);\n\njulia> dmc_est = fit(DiscreteMarkovChain, states);\n\njulia> round.(dmc_est.P, digits=1)\n2×2 Matrix{Float64}:\n 0.9  0.1\n 0.2  0.8","category":"page"},{"location":"tutorial/#Continuous-time-Markov-chains","page":"Tutorial","title":"Continuous time Markov chains","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> cmc = ContinuousMarkovChain(π0 = [0.3, 0.7], Q = [-1. 1.; 2. -2.]);\n\njulia> history = rand(cmc, 0., 1000.);\n\njulia> cmc_est = fit(ContinuousMarkovChain, history);\n\njulia> round.(cmc_est.Q, digits=1)\n2×2 Matrix{Float64}:\n -1.0   1.0\n  1.9  -1.9","category":"page"},{"location":"tutorial/#Hidden-Markov-models","page":"Tutorial","title":"Hidden Markov models","text":"","category":"section"},{"location":"tutorial/#Working-with-point-processes","page":"Tutorial","title":"Working with point processes","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We finally demonstrate the main goal of the package: point process simulation and inference. All point processes are subtypes of AbstractPointProcess{L,M}, where L is the type of event locations and M is the type of event marks.","category":"page"},{"location":"tutorial/#Temporal-Poisson-processes","page":"Tutorial","title":"Temporal Poisson processes","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We provide a number of built-in models, including the basic Poisson process on the real line.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> pp = TemporalPoissonProcess(λ = [0.5, 1., 2.]);\n\njulia> history = rand(pp, 0., 1000.);\n\njulia> pp_init = TemporalPoissonProcess(λ = [1., 1., 1.]);\n\njulia> pp_est = fit(pp_init, history);\n\njulia> round.(pp_est.λ, digits=1)\n3-element Vector{Float64}:\n 0.5\n 1.0\n 2.0","category":"page"},{"location":"tutorial/#Implementing-your-own-temporal-point-process","page":"Tutorial","title":"Implementing your own temporal point process","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To implement your own process, you only have to define a subtype of TemporalPointProcess and write the necessary methods: intensity, mark_distribution, ground_intensity and ground_intensity_bound.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"If these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations.","category":"page"},{"location":"utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utils/","page":"Utilities","title":"Utilities","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"utils/plot.jl\", \"utils/utils.jl\"]","category":"page"}]
}
