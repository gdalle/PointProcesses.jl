var documenterSearchIndex = {"docs":
[{"location":"tuto/markov/#Markov-chains","page":"Markov chains","title":"Markov chains","text":"","category":"section"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"DocTestSetup = quote\n    using PointProcesses\nend","category":"page"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"julia> using Random\n\njulia> Random.seed!(63);","category":"page"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"Some point processes are based on underlying Markov processes, which is why we provide a basic implementation for them.","category":"page"},{"location":"tuto/markov/#Discrete-time","page":"Markov chains","title":"Discrete time","text":"","category":"section"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"julia> dmc = DiscreteMarkovChain(π0 = [0.3, 0.7], P = [0.9 0.1; 0.2 0.8]);\n\njulia> states = rand(dmc, 1000);\n\njulia> dmc_est = fit(DiscreteMarkovChain, states);\n\njulia> round.(dmc_est.P, digits=2)\n2×2 Matrix{Float64}:\n 0.9   0.1\n 0.25  0.75","category":"page"},{"location":"tuto/markov/#Continuous-time","page":"Markov chains","title":"Continuous time","text":"","category":"section"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"julia> cmc = ContinuousMarkovChain(π0 = [0.3, 0.7], Q = [-1. 1.; 2. -2.]);\n\njulia> history = rand(cmc, 0., 1000.);\n\njulia> cmc_est = fit(ContinuousMarkovChain, history);\n\njulia> round.(cmc_est.Q, digits=2)\n2×2 Matrix{Float64}:\n -1.0    1.0\n  1.92  -1.92","category":"page"},{"location":"tuto/models/#Built-in-models","page":"Built-in models","title":"Built-in models","text":"","category":"section"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"DocTestSetup = quote\n    using PointProcesses\nend","category":"page"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"julia> using Random\n\njulia> Random.seed!(63);","category":"page"},{"location":"tuto/models/#Poisson-processes","page":"Built-in models","title":"Poisson processes","text":"","category":"section"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"We provide a number of built-in models, starting with basic Poisson processes on the real line.","category":"page"},{"location":"tuto/models/#Multivariate-Poisson-processes","page":"Built-in models","title":"Multivariate Poisson processes","text":"","category":"section"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"julia> pp = MultivariatePoissonProcess(λ = [0.5, 1., 2.]);\n\njulia> history = rand(pp, 0., 1000.);\n\njulia> pp_est = fit(MultivariatePoissonProcess, history);\n\njulia> round.(pp_est.λ, digits=2)\n3-element Vector{Float64}:\n 0.49\n 0.96\n 2.0","category":"page"},{"location":"tuto/models/#General-Poisson-processes","page":"Built-in models","title":"General Poisson processes","text":"","category":"section"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"Todo","category":"page"},{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/#Event-history","page":"API reference","title":"Event history","text":"","category":"section"},{"location":"api/#Abstract-type","page":"API reference","title":"Abstract type","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"history/abstract.jl\"]","category":"page"},{"location":"api/#PointProcesses.AbstractHistory","page":"API reference","title":"PointProcesses.AbstractHistory","text":"AbstractHistory{L, M}\n\nAbstract supertype for event histories with locations of type L and marks of type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#Temporal-history","page":"API reference","title":"Temporal history","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"history/temporal.jl\"]","category":"page"},{"location":"api/#PointProcesses.TemporalHistory","page":"API reference","title":"PointProcesses.TemporalHistory","text":"TemporalHistory{M}\n\nLinear event histories with marks of type M.\n\nFields\n\ntimes::Vector{Float64}: vector of event times\nmarks::Vector{M}: vector of event marks\ntmin::Float64: start time\ntmax::Float64: end time\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.append!-Tuple{TemporalHistory, TemporalHistory}","page":"API reference","title":"Base.append!","text":"append!(h1::TemporalHistory, h2::TemporalHistory)\n\nAdd all the events of h2 at the end of h1.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.push!-Tuple{TemporalHistory, Any, Any}","page":"API reference","title":"Base.push!","text":"push!(h::TemporalHistory{M}, t::Float64, m::M)\n\nAdd event (t, m) at the end of history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.duration-Tuple{TemporalHistory}","page":"API reference","title":"PointProcesses.duration","text":"duration(h::TemporalHistory)\n\nCompute the difference h.tmax - h.tmin.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.has_events","page":"API reference","title":"PointProcesses.has_events","text":"has_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)\n\nCheck the presence of events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"function"},{"location":"api/#PointProcesses.nb_events","page":"API reference","title":"PointProcesses.nb_events","text":"nb_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)\n\nCount events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"function"},{"location":"api/#PointProcesses.time_change-Tuple{TemporalHistory, Any}","page":"API reference","title":"PointProcesses.time_change","text":"time_change(h, Λ)\n\nApply the time rescaling t mapsto Lambda(t) to history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#General-framework-for-point-processes","page":"API reference","title":"General framework for point processes","text":"","category":"section"},{"location":"api/#Abstract-type-2","page":"API reference","title":"Abstract type","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"point_processes/abstract.jl\"]","category":"page"},{"location":"api/#PointProcesses.AbstractPointProcess","page":"API reference","title":"PointProcesses.AbstractPointProcess","text":"AbstractPointProcess{L, M}\n\nThe common supertype of all point processes with location type L and mark type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.build_transform","page":"API reference","title":"PointProcesses.build_transform","text":"build_transform(pp)\n\nReturn a transformation object from TransformVariables that can turn a Vector{Float64} into a NamedTuple with fields matching those of pp.\n\n\n\n\n\n","category":"function"},{"location":"api/#Temporal-point-processes","page":"API reference","title":"Temporal point processes","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"point_processes/temporal.jl\"]","category":"page"},{"location":"api/#PointProcesses.BoundedTemporalPointProcess","page":"API reference","title":"PointProcesses.BoundedTemporalPointProcess","text":"BoundedTemporalPointProcess{L<:Real, M}\n\nStore a temporal point process P with pre-defined start and end times.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.TemporalPointProcess","page":"API reference","title":"PointProcesses.TemporalPointProcess","text":"TemporalPointProcess{L<:Real, M}\n\nThe common supertype of all temporal point processes (i.e. point processes on the real line) with mark type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.rand-Union{Tuple{M}, Tuple{Random.AbstractRNG, TemporalPointProcess{M}, Any, Any}} where M","page":"API reference","title":"Base.rand","text":"rand(rng, pp, tmin, tmax)\n\nSimulate a temporal point process pp on interval [tmin, tmax) using Ogata's algorithm[Ogata_1981].\n\n[Ogata_1981]: Ogata, Y. (1981), “On Lewis’ simulation method for point processes,” IEEE Transactions on Information Theory, 27, 23–31. https://doi.org/10.1109/TIT.1981.1056305.\n\n\n\n\n\n","category":"method"},{"location":"api/#Distributions.logpdf-Tuple{TemporalPointProcess, TemporalHistory}","page":"API reference","title":"Distributions.logpdf","text":"logpdf(pp, h)\n\nCompute the log probability density function for a temporal point process pp applied to history h:\n\n    log f(h) = sum_i log lambda(t_i  h_i) - Lambda(h)\n\nThe default method uses a loop over events combined with integrated_ground_intensity, but it should be reimplemented for specific processes if faster computation is possible.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.check_residuals-Tuple{TemporalPointProcess, TemporalHistory}","page":"API reference","title":"PointProcesses.check_residuals","text":"check_residuals(pp, h)\n\nCheck whether the point process pp is a good fit for history h by applying Ogata's time rescaling method: if (t_i)_i is a temporal point process with intensity lambda(t), then (Lambda(t_i))_i is a standard temporal Poisson process.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.ground_intensity-Tuple{TemporalPointProcess, TemporalHistory, Any}","page":"API reference","title":"PointProcesses.ground_intensity","text":"ground_intensity(pp, h, t)\n\nCompute the ground intensity for a temporal point process pp applied to history h at time t.\n\nThe ground intensity quantifies the instantaneous risk of an event with any mark occurring at time t[Rasmussen_2018]. It can be expressed as\n\n    lambda_g(th) = sum_m in mathcal lambda(t mh)\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.ground_intensity_bound-Tuple{TemporalPointProcess, TemporalHistory, Any}","page":"API reference","title":"PointProcesses.ground_intensity_bound","text":"ground_intensity_bound(pp, θ, h, t)\n\nCompute a local upper bound on the ground intensity for a temporal point process pp applied to history h at time t[Rasmussen_2018].\n\nReturn a tuple of the form (B L) satisfying\n\n    forall u in t t+L) quad lambda_g(th) leq B\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.integrated_ground_intensity","page":"API reference","title":"PointProcesses.integrated_ground_intensity","text":"integrated_ground_intensity(pp, h[, t])\n\nCompute the integrated ground intensity (or compensator) Lambda(th) for a temporal point process pp applied to history h:\n\n    Lambda(h) = int lambda_g(th) mathrmdt\n\nThe default method uses Quadrature.jl for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.\n\n\n\n\n\n","category":"function"},{"location":"api/#PointProcesses.intensity-Tuple{TemporalPointProcess, TemporalHistory, Any, Any}","page":"API reference","title":"PointProcesses.intensity","text":"intensity(pp, h, t, m)\n\nCompute the conditional intensity for a temporal point process pp applied to history h and event (t, m).\n\nThe conditional intensity function lambda(t m  h) quantifies the instantaneous risk  of an event with mark m occurring at time t[Rasmussen_2018].\n\n[Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.mark_distribution-Tuple{TemporalPointProcess, TemporalHistory, Any}","page":"API reference","title":"PointProcesses.mark_distribution","text":"mark_distribution(pp, h, t)\n\nCompute the distribution of marks for a temporal point process pp knowing that an event takes place at time t after history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#StatsBase.fit-Union{Tuple{PP}, Tuple{M}, Tuple{PP, TemporalHistory{M}}} where {M, PP<:TemporalPointProcess}","page":"API reference","title":"StatsBase.fit","text":"fit(pp0, h)\n\nCompute the optimal parameter for a temporal point process of type typeof(pp0) on history h using maximum likelihood:\n\n    hattheta = mathrmargmax f_theta(h) theta in Theta\n\nThe default method uses GalacticOptim.jl for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.\n\n\n\n\n\n","category":"method"},{"location":"api/#Built-in-models","page":"API reference","title":"Built-in models","text":"","category":"section"},{"location":"api/#Poisson-processes","page":"API reference","title":"Poisson processes","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"models/poisson.jl\", \"models/poisson_inhomogeneous.jl\"]","category":"page"},{"location":"api/#PointProcesses.PoissonProcess","page":"API reference","title":"PointProcesses.PoissonProcess","text":"PoissonProcess{R<:Real,D,M}\n\nHomogeneous temporal Poisson process with arbitrary mark distribution.\n\nFields\n\nλ::R: event rate.\nmark_dist: mark distribution.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.InhomogeneousPoissonProcess","page":"API reference","title":"PointProcesses.InhomogeneousPoissonProcess","text":"InhomogeneousPoissonProcess{D,M}\n\nInhomogeneous temporal Poisson process with arbitrary mark distribution.\n\nFields\n\nλ::Function: intensity function.\nmark_dist::D: mark distribution.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"models/poisson_multivariate.jl\", \"models/poisson_multivariate_naive.jl\"]","category":"page"},{"location":"api/#PointProcesses.MultivariatePoissonProcess","page":"API reference","title":"PointProcesses.MultivariatePoissonProcess","text":"MultivariatePoissonProcess{R}\n\nHomogeneous temporal multivariate Poisson process.\n\nFields\n\nλ::Vector{R}: event rates.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.NaiveMultivariatePoissonProcess","page":"API reference","title":"PointProcesses.NaiveMultivariatePoissonProcess","text":"NaiveMultivariatePoissonProcess{R}\n\nHomogeneous temporal multivariate Poisson process.\n\nThis naive implementation demonstrates the use of our general TemporalPointProcess interface.\n\nFields\n\nλ::Vector{R}: event rates.\n\n\n\n\n\n","category":"type"},{"location":"api/#Hawkes-processes","page":"API reference","title":"Hawkes processes","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"models/hawkes.jl\"]","category":"page"},{"location":"api/#PointProcesses.MultivariateHawkesProcess","page":"API reference","title":"PointProcesses.MultivariateHawkesProcess","text":"MultivariateHawkesProcess{R}\n\nHomogeneous temporal multivariate Hawkes process.\n\nFields\n\nλ::Vector{R}: base event rates\nα::Matrix{R}: excitation amplitudes\nβ::Matrix{R}: excitation decays\n\n\n\n\n\n","category":"type"},{"location":"api/#Markov-chains","page":"API reference","title":"Markov chains","text":"","category":"section"},{"location":"api/#Abstract-type-3","page":"API reference","title":"Abstract type","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"markov/abstract.jl\"]","category":"page"},{"location":"api/#Discrete-time","page":"API reference","title":"Discrete time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"markov/discrete_time.jl\"]","category":"page"},{"location":"api/#PointProcesses.DiscreteMarkovChain","page":"API reference","title":"PointProcesses.DiscreteMarkovChain","text":"DiscreteMarkovChain\n\nDiscrete-time Markov chain with finite state space.\n\nFields\n\nπ0: initial state distribution\nP: state transition matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.DiscreteMarkovChainStats","page":"API reference","title":"PointProcesses.DiscreteMarkovChainStats","text":"DiscreteMarkovChainStats\n\nSufficient statistics for the likelihood of a DiscreteMarkovChain.\n\n\n\n\n\n","category":"type"},{"location":"api/#Continuous-time","page":"API reference","title":"Continuous time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"markov/continuous_time.jl\"]","category":"page"},{"location":"api/#PointProcesses.ContinuousMarkovChain","page":"API reference","title":"PointProcesses.ContinuousMarkovChain","text":"ContinuousMarkovChain\n\nContinuous-time Markov chain with finite state space.\n\nFields\n\nπ0: initial state distribution\nQ: rate matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.ContinuousMarkovChainStats","page":"API reference","title":"PointProcesses.ContinuousMarkovChainStats","text":"ContinuousMarkovChainStats\n\nSufficient statistics for the likelihood of a ContinuousMarkovChain.\n\n\n\n\n\n","category":"type"},{"location":"api/#Hidden-Markov-models","page":"API reference","title":"Hidden Markov models","text":"","category":"section"},{"location":"api/#Discrete-time-2","page":"API reference","title":"Discrete time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"hmm/hmm.jl\", \"hmm/baum_welch.jl\"]","category":"page"},{"location":"api/#PointProcesses.HiddenMarkovModel","page":"API reference","title":"PointProcesses.HiddenMarkovModel","text":"HiddenMarkovModel{Tr<:DiscreteMarkovChain,Em}\n\nHidden Markov Model with arbitrary transition model (must be a discrete Markov chain) and emission distributions.\n\nFields\n\ntransitions::Tr: state evolution process.\nemissions::Vector{Em}: one emission distribution per state.\n\n\n\n\n\n","category":"type"},{"location":"api/#Continuous-time-2","page":"API reference","title":"Continuous time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"hmm/mmpp.jl\", \"hmm/ryden.jl\"]","category":"page"},{"location":"api/#PointProcesses.MarkovModulatedPoissonProcess","page":"API reference","title":"PointProcesses.MarkovModulatedPoissonProcess","text":"MarkovModulatedPoissonProcess{M,Tr<:ContinuousMarkovChain,Em<:PoissonProcess{M}}\n\nMarkov-Modulated Poisson Process with mark type M.\n\nFields\n\ntransitions::Tr: state evolution process.\nemissions::Vector{Em}: one emission distribution per state.\n\n\n\n\n\n","category":"type"},{"location":"api/#Utilities","page":"API reference","title":"Utilities","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"utils/plot.jl\", \"utils/utils.jl\"]","category":"page"},{"location":"tuto/utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"tuto/utils/","page":"Utilities","title":"Utilities","text":"Todo","category":"page"},{"location":"tuto/hmm/#Hidden-Markov-models","page":"Hidden Markov models","title":"Hidden Markov models","text":"","category":"section"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"DocTestSetup = quote\n    using PointProcesses\nend","category":"page"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"julia> using Random\n\njulia> using Distributions\n\njulia> Random.seed!(63);","category":"page"},{"location":"tuto/hmm/#Discrete-time","page":"Hidden Markov models","title":"Discrete time","text":"","category":"section"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"Here is an example of the Baum-Welch estimation algorithm applied to a discrete HMM.","category":"page"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"julia> hmm = HiddenMarkovModel(\n           transitions = DiscreteMarkovChain(π0 = [0.3, 0.7], P = [0.9 0.1; 0.2 0.8]),\n           emissions = [Normal(1, 0.3), Normal(2, 0.3)]\n       );\n\njulia> states, observations = rand(hmm, 1000);\n\njulia> hmm_init = HiddenMarkovModel(\n           transitions = DiscreteMarkovChain(π0 = randprobvec(2), P = randtransmat(2)),\n           emissions = [Normal(rand(), 1), Normal(rand(), 1)]\n       );\n\njulia> hmm_est, logL_evolution = baum_welch(hmm_init, observations, iterations=100);\n\njulia> minimum(diff(logL_evolution)) > -1e-10\ntrue\n\njulia> round.(transition_matrix(hmm_est), digits=2)\n2×2 Matrix{Float64}:\n 0.9   0.1\n 0.23  0.77\n\njulia> round(mean(emission(hmm_est, 1)), digits=2)\n1.02\n\njulia> round(mean(emission(hmm_est, 2)), digits=2)\n1.99","category":"page"},{"location":"roadmap/#Roadmap","page":"Roadmap","title":"Roadmap","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Here is a list of features that are yet to be implemented...","category":"page"},{"location":"roadmap/#Soon","page":"Roadmap","title":"Soon","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Compatibility with reverse-mode AD (Zygote.jl)\nHawkes processes\nCox processes\nMarkov-Modulated Poisson Processes\nPrediction and evaluation utilities","category":"page"},{"location":"roadmap/#Someday","page":"Roadmap","title":"Someday","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Poisson processes on arbitrary measured spaces\nPiecewise-Deterministic Markov Processes\nGaussian Process-modulated Poisson Processes","category":"page"},{"location":"roadmap/#Maybe-not","page":"Roadmap","title":"Maybe not","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Neural point processes?\nInterface with MLJ.jl?\nInterface with Turing.jl?","category":"page"},{"location":"tuto/point_processes/#General-framework-for-point-processes","page":"General framework","title":"General framework for point processes","text":"","category":"section"},{"location":"tuto/point_processes/","page":"General framework","title":"General framework","text":"DocTestSetup = quote\n    using PointProcesses\nend","category":"page"},{"location":"tuto/point_processes/#Working-with-temporal-point-processes","page":"General framework","title":"Working with temporal point processes","text":"","category":"section"},{"location":"tuto/point_processes/","page":"General framework","title":"General framework","text":"The main goal of this package lies in point process simulation and inference. All point processes are subtypes of AbstractPointProcess{L,M}, where L is the type of event locations and M is the type of event marks.","category":"page"},{"location":"tuto/point_processes/","page":"General framework","title":"General framework","text":"To implement your own process, you only have to define a subtype of TemporalPointProcess and write the necessary methods: intensity, mark_distribution, ground_intensity and ground_intensity_bound.","category":"page"},{"location":"tuto/point_processes/","page":"General framework","title":"General framework","text":"As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations.","category":"page"},{"location":"tuto/point_processes/","page":"General framework","title":"General framework","text":"As an example, we included a naive re-implementation of  MultivariatePoissonProcess, called NaiveMultivariatePoissonProcess. Looking at its source may help you understand the requirements of the interface.","category":"page"},{"location":"#PointProcesses","page":"Home","title":"PointProcesses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation of PointProcesses.","category":"page"},{"location":"","page":"Home","title":"Home","text":"PointProcesses","category":"page"},{"location":"#PointProcesses","page":"Home","title":"PointProcesses","text":"A package for point process modeling, simulation and inference.\n\n\n\n\n\n","category":"module"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To get started with code examples, take a look at the Tutorial section. The objects and methods defined in this package are documented in the API reference.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"tuto/history/#Event-history","page":"Event history","title":"Event history","text":"","category":"section"},{"location":"tuto/history/","page":"Event history","title":"Event history","text":"DocTestSetup = quote\n    using PointProcesses\nend","category":"page"},{"location":"tuto/history/","page":"Event history","title":"Event history","text":"To analyze point processes, we need a way to store their realizations. This is the purpose of the AbstractHistory subtypes, of which only TemporalHistory is implemented for now.","category":"page"},{"location":"tuto/history/","page":"Event history","title":"Event history","text":"julia> history = TemporalHistory([0.2, 0.8, 1.1], [\"a\", \"b\", \"c\"], 0.0, 2.0);\n\njulia> duration(history)\n2.0\n\njulia> nb_events(history)\n3\n\njulia> nb_events(history, 1.0, 2.0)\n1\n\njulia> has_events(history)\ntrue\n\njulia> has_events(history, 1.5, 2.0)\nfalse\n\njulia> push!(history, 1.7, \"d\")\n\njulia> has_events(history, 1.5, 2.0)\ntrue","category":"page"}]
}
