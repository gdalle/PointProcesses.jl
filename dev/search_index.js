var documenterSearchIndex = {"docs":
[{"location":"tuto/markov/#Markov-chains","page":"Markov chains","title":"Markov chains","text":"","category":"section"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"Some point processes are based on underlying Markov processes, which is why we provide a basic implementation for them.","category":"page"},{"location":"tuto/markov/#Discrete-time","page":"Markov chains","title":"Discrete time","text":"","category":"section"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"using PointProcesses, Random; Random.seed!(63)\n\ndmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])\nstates = rand(dmc, 1000);\ndmc_est = fit(DiscreteMarkovChain, states)\ntransition_matrix(dmc_est)","category":"page"},{"location":"tuto/markov/#Continuous-time","page":"Markov chains","title":"Continuous time","text":"","category":"section"},{"location":"tuto/markov/","page":"Markov chains","title":"Markov chains","text":"using PointProcesses, Random; Random.seed!(63)\n\ncmc = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])\nhistory = rand(cmc, 0., 1000.);\ncmc_est = fit(ContinuousMarkovChain, history)\nrate_matrix(cmc_est)","category":"page"},{"location":"tuto/models/#Built-in-models","page":"Built-in models","title":"Built-in models","text":"","category":"section"},{"location":"tuto/models/#Poisson-processes","page":"Built-in models","title":"Poisson processes","text":"","category":"section"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"We provide a general PoissonProcess structure with a single intensity value and an arbitrary mark distribution.","category":"page"},{"location":"tuto/models/","page":"Built-in models","title":"Built-in models","text":"using MeasureTheory, PointProcesses, Random; Random.seed!(63)\n\nmark_dist = Dists.Categorical([0.1, 0.3, 0.6])\npp = PoissonProcess(5., mark_dist)\nh = rand(pp, 0., 100.)\npp_est = fit(PoissonProcess{Dists.Categorical}, h)\nintensity(pp_est)\nDists.probs(mark_distribution(pp_est))","category":"page"},{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/#Event-history","page":"API reference","title":"Event history","text":"","category":"section"},{"location":"api/#Abstract-type","page":"API reference","title":"Abstract type","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"history/abstract.jl\"]","category":"page"},{"location":"api/#PointProcesses.AbstractHistory","page":"API reference","title":"PointProcesses.AbstractHistory","text":"AbstractHistory{L,M}\n\nAbstract supertype for event histories with locations of type L and marks of type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#Temporal-history","page":"API reference","title":"Temporal history","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"history/temporal.jl\"]","category":"page"},{"location":"api/#PointProcesses.TemporalHistory","page":"API reference","title":"PointProcesses.TemporalHistory","text":"TemporalHistory{M}\n\nLinear event histories with marks of type M.\n\nFields\n\ntimes::Vector{Float64}: vector of event times\nmarks::Vector{M}: vector of event marks\ntmin::Float64: start time\ntmax::Float64: end time\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.append!-Tuple{TemporalHistory, TemporalHistory}","page":"API reference","title":"Base.append!","text":"append!(h1::TemporalHistory, h2::TemporalHistory)\n\nAdd all the events of h2 at the end of h1.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.push!-Tuple{TemporalHistory, Any, Any}","page":"API reference","title":"Base.push!","text":"push!(h::TemporalHistory{M}, t::Float64, m::M)\n\nAdd event (t, m) at the end of history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.duration-Tuple{TemporalHistory}","page":"API reference","title":"PointProcesses.duration","text":"duration(h::TemporalHistory)\n\nCompute the difference h.tmax - h.tmin.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.has_events","page":"API reference","title":"PointProcesses.has_events","text":"has_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)\n\nCheck the presence of events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"function"},{"location":"api/#PointProcesses.max_time-Tuple{TemporalHistory}","page":"API reference","title":"PointProcesses.max_time","text":"max_time(h)\n\nReturn the end time of h (not the same as the last event time).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.min_time-Tuple{TemporalHistory}","page":"API reference","title":"PointProcesses.min_time","text":"min_time(h)\n\nReturn the starting time of h (not the same as the first event time).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.nb_events","page":"API reference","title":"PointProcesses.nb_events","text":"nb_events(h::TemporalHistory, tmin=-Inf, tmax=Inf)\n\nCount events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"function"},{"location":"api/#PointProcesses.time_change-Tuple{TemporalHistory, Any}","page":"API reference","title":"PointProcesses.time_change","text":"time_change(h, Λ)\n\nApply the time rescaling t mapsto Lambda(t) to history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#Point-processes","page":"API reference","title":"Point processes","text":"","category":"section"},{"location":"api/#Abstract-type-2","page":"API reference","title":"Abstract type","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"point_processes/abstract.jl\"]","category":"page"},{"location":"api/#PointProcesses.AbstractPointProcess","page":"API reference","title":"PointProcesses.AbstractPointProcess","text":"AbstractPointProcess{L, M}\n\nThe common supertype of all point processes with location type L and mark type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.build_transform","page":"API reference","title":"PointProcesses.build_transform","text":"build_transform(pp)\n\nReturn a transformation object from TransformVariables that can turn a Vector{Float64} into a NamedTuple with fields matching those of pp.\n\n\n\n\n\n","category":"function"},{"location":"api/#Temporal-point-processes","page":"API reference","title":"Temporal point processes","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"point_processes/temporal.jl\"]","category":"page"},{"location":"api/#PointProcesses.BoundedTemporalPointProcess","page":"API reference","title":"PointProcesses.BoundedTemporalPointProcess","text":"BoundedTemporalPointProcess{L<:Real, M}\n\nStore a temporal point process P with pre-defined start and end times.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.TemporalPointProcess","page":"API reference","title":"PointProcesses.TemporalPointProcess","text":"TemporalPointProcess{L<:Real, M}\n\nThe common supertype of all temporal point processes (i.e. point processes on the real line) with mark type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.rand-Union{Tuple{M}, Tuple{Random.AbstractRNG, TemporalPointProcess{M}, Any, Any}} where M","page":"API reference","title":"Base.rand","text":"rand(rng, pp, tmin, tmax)\n\nSimulate a temporal point process pp on interval [tmin, tmax) using Ogata's algorithm[Ogata_1981].\n\n[Ogata_1981]: Ogata, Y. (1981), “On Lewis’ simulation method for point processes,” IEEE Transactions on Information Theory, 27, 23–31. https://doi.org/10.1109/TIT.1981.1056305.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.check_residuals-Tuple{TemporalPointProcess, TemporalHistory}","page":"API reference","title":"PointProcesses.check_residuals","text":"check_residuals(pp, h)\n\nCheck whether the point process pp is a good fit for history h by applying Ogata's time rescaling method: if (t_i)_i is a temporal point process with intensity lambda(t), then (Lambda(t_i))_i is a standard temporal Poisson process.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.ground_intensity-Tuple{TemporalPointProcess, TemporalHistory, Any}","page":"API reference","title":"PointProcesses.ground_intensity","text":"ground_intensity(pp, h, t)\n\nCompute the ground intensity for a temporal point process pp applied to history h at time t.\n\nThe ground intensity quantifies the instantaneous risk of an event with any mark occurring at time t[Rasmussen_2018]. It can be expressed as\n\n    lambda_g(th) = sum_m in mathcal lambda(t mh)\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.ground_intensity_bound-Tuple{TemporalPointProcess, TemporalHistory, Any}","page":"API reference","title":"PointProcesses.ground_intensity_bound","text":"ground_intensity_bound(pp, θ, h, t)\n\nCompute a local upper bound on the ground intensity for a temporal point process pp applied to history h at time t[Rasmussen_2018].\n\nReturn a tuple of the form (B L) satisfying\n\n    forall u in t t+L) quad lambda_g(th) leq B\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.integrated_ground_intensity","page":"API reference","title":"PointProcesses.integrated_ground_intensity","text":"integrated_ground_intensity(pp, h[, t])\n\nCompute the integrated ground intensity (or compensator) Lambda(th) for a temporal point process pp applied to history h:\n\n    Lambda(h) = int lambda_g(th) mathrmdt\n\nThe default method uses Quadrature.jl for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.\n\n\n\n\n\n","category":"function"},{"location":"api/#PointProcesses.intensity-Tuple{TemporalPointProcess, TemporalHistory, Any, Any}","page":"API reference","title":"PointProcesses.intensity","text":"intensity(pp, h, t, m)\n\nCompute the conditional intensity for a temporal point process pp applied to history h and event (t, m).\n\nThe conditional intensity function lambda(t m  h) quantifies the instantaneous risk  of an event with mark m occurring at time t[Rasmussen_2018].\n\n[Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.logpdf-Tuple{TemporalPointProcess, TemporalHistory}","page":"API reference","title":"PointProcesses.logpdf","text":"logpdf(pp, h)\n\nCompute the log probability density function for a temporal point process pp applied to history h:\n\n    log f(h) = sum_i log lambda(t_i  h_i) - Lambda(h)\n\nThe default method uses a loop over events combined with integrated_ground_intensity, but it should be reimplemented for specific processes if faster computation is possible.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.mark_distribution-Tuple{TemporalPointProcess, TemporalHistory, Any}","page":"API reference","title":"PointProcesses.mark_distribution","text":"mark_distribution(pp, h, t)\n\nCompute the distribution of marks for a temporal point process pp knowing that an event takes place at time t after history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#StatsBase.fit-Union{Tuple{PP}, Tuple{M}, Tuple{PP, TemporalHistory{M}}} where {M, PP<:TemporalPointProcess}","page":"API reference","title":"StatsBase.fit","text":"fit(pp0, h)\n\nCompute the optimal parameter for a temporal point process of type typeof(pp0) on history h using maximum likelihood:\n\n    hattheta = mathrmargmax f_theta(h) theta in Theta\n\nThe default method uses GalacticOptim.jl for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.\n\n\n\n\n\n","category":"method"},{"location":"api/#Built-in-models","page":"API reference","title":"Built-in models","text":"","category":"section"},{"location":"api/#Poisson-processes","page":"API reference","title":"Poisson processes","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"models/poisson.jl\", \"models/poisson_multivariate_naive.jl\"]","category":"page"},{"location":"api/#PointProcesses.PoissonProcess","page":"API reference","title":"PointProcesses.PoissonProcess","text":"PoissonProcess{D,M}\n\nHomogeneous temporal Poisson process with arbitrary mark distribution.\n\nFields\n\nλ::Float64: event rate.\nmark_dist::D: mark distribution.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.NaiveMultivariatePoissonProcess","page":"API reference","title":"PointProcesses.NaiveMultivariatePoissonProcess","text":"NaiveMultivariatePoissonProcess{R}\n\nHomogeneous temporal multivariate Poisson process.\n\nThis naive implementation demonstrates the use of our general TemporalPointProcess interface.\n\nFields\n\nλ::Vector{R}: event rates.\n\n\n\n\n\n","category":"type"},{"location":"api/#Hawkes-processes","page":"API reference","title":"Hawkes processes","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"models/hawkes.jl\"]","category":"page"},{"location":"api/#Markov-chains","page":"API reference","title":"Markov chains","text":"","category":"section"},{"location":"api/#Abstract-type-3","page":"API reference","title":"Abstract type","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"markov/abstract.jl\"]","category":"page"},{"location":"api/#Discrete-time","page":"API reference","title":"Discrete time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"markov/discrete_time.jl\"]","category":"page"},{"location":"api/#PointProcesses.DiscreteMarkovChain","page":"API reference","title":"PointProcesses.DiscreteMarkovChain","text":"DiscreteMarkovChain{T1,T2}\n\nDiscrete-time Markov chain with finite state space.\n\nFields\n\nπ0::T1: initial state distribution.\nP::T2: state transition matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.DiscreteMarkovChainPrior","page":"API reference","title":"PointProcesses.DiscreteMarkovChainPrior","text":"DiscreteMarkovChainPrior\n\nDefine a Dirichlet prior on the initial distribution and on the transitions from each state.\n\nFields\n\nπ0α: Dirichlet parameter for the initial distribution\nPα: Dirichlet parameters for the transition matrix\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.DiscreteMarkovChainStats","page":"API reference","title":"PointProcesses.DiscreteMarkovChainStats","text":"DiscreteMarkovChainStats\n\nSufficient statistics for the likelihood of a DiscreteMarkovChain.\n\n\n\n\n\n","category":"type"},{"location":"api/#Continuous-time","page":"API reference","title":"Continuous time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"markov/continuous_time.jl\"]","category":"page"},{"location":"api/#PointProcesses.ContinuousMarkovChain","page":"API reference","title":"PointProcesses.ContinuousMarkovChain","text":"ContinuousMarkovChain{T1,T2}\n\nContinuous-time Markov chain with finite state space.\n\nFields\n\nπ0::T1: initial state distribution\nQ::T2: rate matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.ContinuousMarkovChainPrior","page":"API reference","title":"PointProcesses.ContinuousMarkovChainPrior","text":"ContinuousMarkovChainPrior\n\nDefine a Dirichlet prior on the initial distribution and a Gamma prior on the transition rates from each state.\n\nFields\n\nπ0α: Dirichlet parameter\nPα: Gamma rate parameters\nPβ: Gamma shape parameters\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.ContinuousMarkovChainStats","page":"API reference","title":"PointProcesses.ContinuousMarkovChainStats","text":"ContinuousMarkovChainStats\n\nSufficient statistics for the likelihood of a ContinuousMarkovChain.\n\n\n\n\n\n","category":"type"},{"location":"api/#Hidden-Markov-models","page":"API reference","title":"Hidden Markov models","text":"","category":"section"},{"location":"api/#Discrete-time-2","page":"API reference","title":"Discrete time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"hmm/hmm.jl\", \"hmm/baum_welch.jl\"]","category":"page"},{"location":"api/#PointProcesses.HiddenMarkovModel","page":"API reference","title":"PointProcesses.HiddenMarkovModel","text":"HMM{Tr<:DiscreteMarkovChain,Em}\n\nHidden Markov Model with arbitrary transition model (must be a discrete Markov chain) and emission distributions.\n\nFields\n\ntransitions::Tr: state evolution process.\nemissions::Vector{Em}: one emission distribution per state.\n\n\n\n\n\n","category":"type"},{"location":"api/#Continuous-time-2","page":"API reference","title":"Continuous time","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nPages = [\"hmm/mmpp.jl\", \"hmm/ryden.jl\"]","category":"page"},{"location":"api/#PointProcesses.MarkovModulatedPoissonProcess","page":"API reference","title":"PointProcesses.MarkovModulatedPoissonProcess","text":"MarkovModulatedPoissonProcess{M,Tr<:ContinuousMarkovChain,Em<:TemporalPointProcess{M}}\n\nMarkov-Modulated Poisson Process with mark type M.\n\nFields\n\ntransitions::Tr: state evolution process.\nemissions::Vector{Em}: one emission distribution per state.\n\n\n\n\n\n","category":"type"},{"location":"api/#Utilities","page":"API reference","title":"Utilities","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]\nOrder = [:type, :function]\nPages = [\"utils/categorical.jl\", \"utils/overflow.jl\", \"utils/plot.jl\", \"utils/randvals.jl\"]","category":"page"},{"location":"api/#PointProcesses.plot_events-Union{Tuple{TemporalHistory{M}}, Tuple{M}} where M<:Real","page":"API reference","title":"PointProcesses.plot_events","text":"plot_events(h)\n\nPlot an event history with its marks if those are of a real subtype.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.plot_intensity-Union{Tuple{M}, Tuple{TemporalPointProcess{M}, TemporalHistory{M}}} where M","page":"API reference","title":"PointProcesses.plot_intensity","text":"plot_intensity(pp, h; npoints)\n\nPlot the conditional intensity function of a temporal point process along a given history.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.qqplot_interevent_times-Tuple{TemporalHistory}","page":"API reference","title":"PointProcesses.qqplot_interevent_times","text":"qqplot_interevent_times(h)\n\nCompare the distribution of inter-event times in a given history to the exponential distribution.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.randprobvec-Tuple{Any}","page":"API reference","title":"PointProcesses.randprobvec","text":"randprobvec(n)\n\nReturn a random probability distribution vector of size n.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.randtransmat-Tuple{Any}","page":"API reference","title":"PointProcesses.randtransmat","text":"randtransmat(n)\n\nReturn a stochastic matrix of size n with random transition probability distributions.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.uniformprobvec-Tuple{Any}","page":"API reference","title":"PointProcesses.uniformprobvec","text":"uniformprobvec(n)\n\nReturn a uniform probability distribution vector of size n.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.uniformtransmat-Tuple{Any}","page":"API reference","title":"PointProcesses.uniformtransmat","text":"uniformtransmat(n)\n\nReturn a stochastic matrix of size n with uniform transition probability distributions.\n\n\n\n\n\n","category":"method"},{"location":"tuto/utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"tuto/utils/","page":"Utilities","title":"Utilities","text":"Todo","category":"page"},{"location":"tuto/hmm/#Hidden-Markov-models","page":"Hidden Markov models","title":"Hidden Markov models","text":"","category":"section"},{"location":"tuto/hmm/#Discrete-time","page":"Hidden Markov models","title":"Discrete time","text":"","category":"section"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"Here is an example of the Baum-Welch estimation algorithm applied to a discrete HMM.","category":"page"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"using MeasureTheory, PointProcesses, Random; Random.seed!(63)\n\ntr = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])\nem1 = Dists.Normal(1, 0.3)\nem2 = Dists.Normal(2, 0.3)\nhmm = HiddenMarkovModel(tr, [em1, em2])\n\nstates, observations = rand(hmm, 1000)\n\ntr_init = DiscreteMarkovChain(randprobvec(2), randtransmat(2))\nem1_init = Dists.Normal(rand(), 1)\nem2_init = Dists.Normal(rand(), 1)\nhmm_init = HiddenMarkovModel(tr_init, [em1_init, em2_init])\n\nhmm_est, logL_evolution = baum_welch(hmm_init, observations, iterations=100)\n\ntransition_matrix(hmm_est)\n\nemissions(hmm_est)","category":"page"},{"location":"tuto/hmm/#Continuous-time","page":"Hidden Markov models","title":"Continuous time","text":"","category":"section"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"Here is an example of the Baum-Welch estimation algorithm applied to a MMPP.","category":"page"},{"location":"tuto/hmm/","page":"Hidden Markov models","title":"Hidden Markov models","text":"using MeasureTheory, PointProcesses, Random; Random.seed!(63)\n\ntr = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])\nem1 = PoissonProcess(1., Dists.Normal(1, 1))\nem2 = PoissonProcess(2., Dists.Normal(-1, 1))\nmmpp = MarkovModulatedPoissonProcess(tr, [em1, em2])\n\nstate_history, observations = rand(mmpp, 0., 1000.)","category":"page"},{"location":"roadmap/#Roadmap","page":"Roadmap","title":"Roadmap","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Here is a list of features that are yet to be implemented...","category":"page"},{"location":"roadmap/#Soon","page":"Roadmap","title":"Soon","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Compatibility with reverse-mode AD (Zygote.jl)\nHawkes processes\nCox processes\nMarkov-Modulated Poisson Processes\nPrediction and evaluation utilities","category":"page"},{"location":"roadmap/#Someday","page":"Roadmap","title":"Someday","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Poisson processes on arbitrary measured spaces\nPiecewise-Deterministic Markov Processes\nGaussian Process-modulated Poisson Processes","category":"page"},{"location":"roadmap/#Maybe-not","page":"Roadmap","title":"Maybe not","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Neural point processes?\nInterface with MLJ.jl?\nInterface with Turing.jl?","category":"page"},{"location":"tuto/point_processes/#General-framework-for-point-processes","page":"Point processes","title":"General framework for point processes","text":"","category":"section"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"using ForwardDiff, GalacticOptim, Optim, Random\nRandom.seed!(63)","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"The main goal of this package lies in point process simulation and inference. All point processes are subtypes of AbstractPointProcess{L,M}, where L is the type of event locations and M is the type of event marks.","category":"page"},{"location":"tuto/point_processes/#Temporal-point-processes","page":"Point processes","title":"Temporal point processes","text":"","category":"section"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"Temporal point processes have specific features which make it possible to define generic functions to simulate or estimate them.[Rasmussen_2018]","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"[Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"To implement your own process, you only have to define a subtype of TemporalPointProcess and implement the necessary methods: intensity, mark_distribution, ground_intensity and ground_intensity_bound. As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations (for instance to avoid numerical optimization or gradient optimization).","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"As an example, we included a naive implementation of a Poisson process with categorical mark distribution, called NaiveMultivariatePoissonProcess. Looking at its source code may help you understand the requirements of the interface. In the meantime, we can check that the implementation works:","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"using GalacticOptim, PointProcesses, Random; Random.seed!(63)\n\npp = NaiveMultivariatePoissonProcess([1., 2., 3.])\nh = rand(pp, 0., 100.)\npp_init = NaiveMultivariatePoissonProcess(ones(3))\npp_est1 = fit(pp_init, h)\nintensity(pp_est1)","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"We can also tune the automatic differentiation method and the optimization algorithm:","category":"page"},{"location":"tuto/point_processes/","page":"Point processes","title":"Point processes","text":"using ForwardDiff, Optim\n\npp_est2 = fit(pp_init, h, adtype=GalacticOptim.AutoForwardDiff(), alg=Optim.LBFGS())\nintensity(pp_est2)","category":"page"},{"location":"#PointProcesses","page":"Home","title":"PointProcesses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation of PointProcesses.","category":"page"},{"location":"","page":"Home","title":"Home","text":"PointProcesses","category":"page"},{"location":"#PointProcesses","page":"Home","title":"PointProcesses","text":"A package for point process modeling, simulation and inference.\n\n\n\n\n\n","category":"module"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the latest release available in the General registry, run the following code in Julia's REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(\"PointProcesses\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"To install the current development version, run this instead:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(\"https://github.com/gdalle/PointProcesses.jl\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to get started with code examples, take a look at the Tutorial section. The objects and functions defined in this package are documented in the API reference.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"tuto/history/#Event-history","page":"Event history","title":"Event history","text":"","category":"section"},{"location":"tuto/history/","page":"Event history","title":"Event history","text":"ENV[\"GKSwstype\"] = \"100\"","category":"page"},{"location":"tuto/history/","page":"Event history","title":"Event history","text":"To analyze point processes, we need a way to store their realizations. This is the purpose of the AbstractHistory subtypes, of which only TemporalHistory is implemented for now.","category":"page"},{"location":"tuto/history/","page":"Event history","title":"Event history","text":"using PointProcesses\n\nhistory = TemporalHistory([0.2, 0.8, 1.1], [\"a\", \"b\", \"c\"], 0.0, 2.0);\nduration(history)\nnb_events(history)\nnb_events(history, 1.0, 2.0)\nhas_events(history)\nhas_events(history, 1.5, 2.0)\npush!(history, 1.7, \"d\")\nhas_events(history, 1.5, 2.0)","category":"page"}]
}
