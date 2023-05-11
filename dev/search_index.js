var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/#Index","page":"API reference","title":"Index","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]","category":"page"},{"location":"api/#Full-docs","page":"API reference","title":"Full docs","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [PointProcesses]","category":"page"},{"location":"api/#PointProcesses.PointProcesses","page":"API reference","title":"PointProcesses.PointProcesses","text":"PointProcesses\n\nA package for temporal point process modeling, simulation and inference.\n\n\n\n\n\n","category":"module"},{"location":"api/#PointProcesses.AbstractPointProcess","page":"API reference","title":"PointProcesses.AbstractPointProcess","text":"AbstractPointProcess{M}\n\nThe common supertype of all temporal point processes (i.e. point processes on the real line) with mark type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.BoundedPointProcess","page":"API reference","title":"PointProcesses.BoundedPointProcess","text":"BoundedPointProcess{M,P,T}\n\nStore a temporal point process P with pre-defined start and end times.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.History","page":"API reference","title":"PointProcesses.History","text":"History{M,T}\n\nLinear event histories with marks of type M and locations of type T, usually times represented as real numbers.\n\nFields\n\ntimes::Vector{T}: sorted vector of event times\nmarks::Vector{M}: associated vector of event marks\ntmin::T: start time\ntmax::T: end time\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.MarkedPoissonProcess","page":"API reference","title":"PointProcesses.MarkedPoissonProcess","text":"MarkedPoissonProcess{R,M,D}\n\nHomogeneous temporal Poisson process with arbitrary mark distribution.\n\nFields\n\nλ::R: event rate.\nmark_dist::D: mark distribution with sample type M.\n\n\n\n\n\n","category":"type"},{"location":"api/#PointProcesses.MultivariatePoissonProcess","page":"API reference","title":"PointProcesses.MultivariatePoissonProcess","text":"MultivariatePoissonProcess{R}\n\nMultivariate homogeneous temporal Poisson process.\n\nFields\n\nλ::Vector{R}: event rates.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.append!-Tuple{History, History}","page":"API reference","title":"Base.append!","text":"append!(h1::History, h2::History)\n\nAdd all the events of h2 at the end of h1.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.push!-Tuple{History, Any, Any}","page":"API reference","title":"Base.push!","text":"push!(h::History, t, m)\n\nAdd event (t, m) at the end of history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#DensityInterface.logdensityof-Tuple{AbstractPointProcess, Any}","page":"API reference","title":"DensityInterface.logdensityof","text":"logdensityof(pp, h)\n\nCompute the log probability density function for a temporal point process pp applied to history h:\n\nℓ(h) = Σₖ log λ(tₖ|hₖ) - Λ(h)\n\nThe default method uses a loop over events combined with integrated_ground_intensity, but it should be reimplemented for specific processes if faster computation is possible.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.check_residuals-Tuple{AbstractPointProcess, Any}","page":"API reference","title":"PointProcesses.check_residuals","text":"check_residuals(pp, h)\n\nCheck whether the point process pp is a good fit for history h by applying Ogata's time rescaling method. If (tₖ)ₖ is a temporal point process with intensity λ, then (Λ(tₖ))ₖ is a standard temporal Poisson process.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.duration-Tuple{History}","page":"API reference","title":"PointProcesses.duration","text":"duration(h::History)\n\nCompute the difference h.tmax - h.tmin.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.event_marks-Tuple{History}","page":"API reference","title":"PointProcesses.event_marks","text":"event_marks(h::History)\n\nReturn the vector of event marks for h, sorted according to their event times.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.event_times-Tuple{History}","page":"API reference","title":"PointProcesses.event_times","text":"event_times(h::History)\n\nReturn the sorted vector of event times for h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.ground_intensity-Union{Tuple{PP}, Tuple{PP, Any, Any}} where PP<:AbstractPointProcess","page":"API reference","title":"PointProcesses.ground_intensity","text":"ground_intensity(pp, h, t)\n\nCompute the ground intensity for a temporal point process pp applied to history h at time t.\n\nThe ground intensity quantifies the instantaneous risk of an event with any mark occurring at time t:\n\nλg(t|h) = Σₘ λ(t,m|h)\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.ground_intensity_bound-Union{Tuple{PP}, Tuple{PP, Any, Any}} where PP<:AbstractPointProcess","page":"API reference","title":"PointProcesses.ground_intensity_bound","text":"ground_intensity_bound(pp, t, h)\n\nCompute a local upper bound on the ground intensity for a temporal point process pp applied to history h at time t.\n\nReturn a tuple of the form (B, L) satisfying λg(t|h) ≤ B for all u ∈ [t, t+L).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.has_events-Tuple{History, Any, Any}","page":"API reference","title":"PointProcesses.has_events","text":"has_events(h::History, tmin, tmax)\n\nCheck the presence of events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.has_events-Tuple{History}","page":"API reference","title":"PointProcesses.has_events","text":"has_events(h::History)\n\nCheck the presence of events in h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.integrated_ground_intensity-Union{Tuple{PP}, Tuple{PP, Any, Any, Any}} where PP<:AbstractPointProcess","page":"API reference","title":"PointProcesses.integrated_ground_intensity","text":"integrated_ground_intensity(pp, h, a, b)\n\nCompute the integrated ground intensity (or compensator) Λ(t|h) for a temporal point process pp applied to history h on interval [a, b]:\n\nΛ(h) = ∫ λg(t|h) dt\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.intensity-Tuple{AbstractPointProcess, Any, Any, Any}","page":"API reference","title":"PointProcesses.intensity","text":"intensity(pp, m, t, h)\n\nCompute the conditional intensity for a temporal point process pp applied to history h and event (t, m).\n\nThe conditional intensity function λ(t,m|h) quantifies the instantaneous risk of an event with mark m occurring at time t after history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.log_intensity-Tuple{AbstractPointProcess, Any, Any, Any}","page":"API reference","title":"PointProcesses.log_intensity","text":"log_intensity(pp, m, t, h)\n\nCompute the logarithm of the conditional intensity for a temporal point process pp applied to history h and event (t, m).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.mark_distribution-Union{Tuple{PP}, Tuple{PP, Any, Any}} where PP<:AbstractPointProcess","page":"API reference","title":"PointProcesses.mark_distribution","text":"mark_distribution(pp, t, h)\n\nCompute the distribution of marks for a temporal point process pp knowing that an event takes place at time t after history h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.max_time-Tuple{History}","page":"API reference","title":"PointProcesses.max_time","text":"max_time(h)\n\nReturn the end time of h (not the same as the last event time).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.min_time-Tuple{History}","page":"API reference","title":"PointProcesses.min_time","text":"min_time(h)\n\nReturn the starting time of h (not the same as the first event time).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.nb_events-Tuple{History}","page":"API reference","title":"PointProcesses.nb_events","text":"nb_events(h::History)\n\nCount events in h.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.nb_events-Union{Tuple{T}, Tuple{M}, Tuple{History{M, T}, Any, Any}} where {M, T}","page":"API reference","title":"PointProcesses.nb_events","text":"nb_events(h::History, tmin, tmax)\n\nCount events in h during the interval [tmin, tmax).\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.simulate_ogata-Union{Tuple{T}, Tuple{M}, Tuple{Random.AbstractRNG, AbstractPointProcess{M}, T, T}} where {M, T<:Real}","page":"API reference","title":"PointProcesses.simulate_ogata","text":"simulate_ogata(rng, pp, tmin, tmax)\n\nSimulate a temporal point process pp on interval [tmin, tmax) using Ogata's algorithm.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.split_into_chunks-Union{Tuple{M}, Tuple{History{M}, Any}} where M","page":"API reference","title":"PointProcesses.split_into_chunks","text":"split_into_chunks(h, chunk_duration)\n\nSplit h into a vector of consecutive histories with individual duration chunk_duration.\n\n\n\n\n\n","category":"method"},{"location":"api/#PointProcesses.time_change-Tuple{History, Any}","page":"API reference","title":"PointProcesses.time_change","text":"time_change(h, Λ)\n\nApply the time rescaling t -> Λ(t) to history h.\n\n\n\n\n\n","category":"method"},{"location":"roadmap/#Roadmap","page":"Roadmap","title":"Roadmap","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Here is a list of features that are yet to be implemented...","category":"page"},{"location":"roadmap/#Soon","page":"Roadmap","title":"Soon","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Compatibility with reverse-mode AD (Zygote.jl)\nHawkes processes\nCox processes\nMarkov-Modulated Poisson Processes\nPrediction and evaluation utilities","category":"page"},{"location":"roadmap/#Someday","page":"Roadmap","title":"Someday","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Poisson processes on arbitrary measured spaces\nPiecewise-Deterministic Markov Processes\nGaussian Process-modulated Poisson Processes","category":"page"},{"location":"roadmap/#Maybe-not","page":"Roadmap","title":"Maybe not","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Neural point processes?\nInterface with MLJ.jl?\nInterface with Turing.jl?","category":"page"},{"location":"#PointProcesses.jl","page":"Home","title":"PointProcesses.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation of PointProcesses.jl, a package for modeling, simulation and inference with temporal point processes.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the latest release available in the General registry, run the following code in Julia's REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.add(\"PointProcesses\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"To install the current development version, run this code instead:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.add(url=\"https://github.com/gdalle/PointProcesses.jl\")","category":"page"}]
}
