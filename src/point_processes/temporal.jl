"""
    TemporalPointProcess{L<:Real, M}

The common supertype of all temporal point processes (i.e. point processes on the real line) with mark type `M`.
"""
abstract type TemporalPointProcess{M} <: AbstractPointProcess{Float64,M} end

"""
    BoundedTemporalPointProcess{L<:Real, M}

Store a temporal point process `P` with pre-defined start and end times.
"""
struct BoundedTemporalPointProcess{M,P<:TemporalPointProcess{M}} <: TemporalPointProcess{M}
    pp::P
    tmin::Float64
    tmax::Float64
end

## Intensity functions

@doc raw"""
    intensity(pp, h, t, m)

Compute the conditional intensity for a temporal point process `pp` applied to history `h` and event `(t, m)`.

The conditional intensity function ``\lambda(t, m | h)`` quantifies the instantaneous risk  of an event with mark `m` occurring at time `t`[^Rasmussen_2018].

[^Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].
"""
function intensity(pp::TemporalPointProcess, h::TemporalHistory, t, m) end

@doc raw"""
    mark_distribution(pp, h, t)

Compute the distribution of marks for a temporal point process `pp` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(pp::TemporalPointProcess, h::TemporalHistory, t) end

@doc raw"""
    ground_intensity(pp, h, t)

Compute the ground intensity for a temporal point process `pp` applied to history `h` at time `t`.

The ground intensity quantifies the instantaneous risk of an event with any mark occurring at time `t`[^Rasmussen_2018]. It can be expressed as
```math
    \lambda_g(t|h) = \sum_{m \in \mathcal} \lambda(t, m|h).
```
"""
function ground_intensity(pp::TemporalPointProcess, h::TemporalHistory, t) end

@doc raw"""
    ground_intensity_bound(pp, θ, h, t)

Compute a local upper bound on the ground intensity for a temporal point process `pp` applied to history `h` at time `t`[^Rasmussen_2018].

Return a tuple of the form $(B, L)$ satisfying
```math
    \forall u \in [t, t+L), \quad \lambda_g(t|h) \leq B
```
"""
function ground_intensity_bound(pp::TemporalPointProcess, h::TemporalHistory, t) end

## Simulation with Ogata's algorithm

Base.eltype(::Type{<:TemporalPointProcess{M}}) where {M} = TemporalHistory{M}

"""
    rand(rng, pp, tmin, tmax)

Simulate a temporal point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm[^Ogata_1981].

[^Ogata_1981]: Ogata, Y. (1981), “On Lewis’ simulation method for point processes,” IEEE Transactions on Information Theory, 27, 23–31. https://doi.org/10.1109/TIT.1981.1056305.
"""
function Base.rand(rng::AbstractRNG, pp::TemporalPointProcess{M}, tmin, tmax) where {M}
    h = TemporalHistory(Float64[], M[], tmin, tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, h, t + eps(t))
        T = B > 0 ? rand(rng, Exponential(1 / B)) : Inf
        if T > L
            t = t + L
        elseif T <= L
            U = rand(Uniform(0, 1))
            if U < ground_intensity(pp, h, t + T) / B
                m = rand(rng, mark_distribution(pp, h, t + T))
                if t + T < tmax
                    push!(h, t + T, m)
                end
            end
            t = t + T
        end
    end
    return h
end

Base.rand(rng::AbstractRNG, tpp::BoundedTemporalPointProcess) =
    rand(rng, tpp.pp, tpp.tmin, tpp.tmax)

Base.rand(tpp::TemporalPointProcess, args...) = rand(Random.GLOBAL_RNG, tpp, args...)

## Learning

@doc raw"""
    integrated_ground_intensity(pp, h[, t])

Compute the integrated ground intensity (or compensator) $\Lambda(t|h)$ for a temporal point process `pp` applied to history `h`:
```math
    \Lambda(h) = \int \lambda_g(t|h) \mathrm{d}t.
```

The default method uses [Quadrature.jl](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.
"""
function integrated_ground_intensity(
    pp::TemporalPointProcess,
    h::TemporalHistory,
    t = max_time(h),
)
    par = [pp]
    f = (t, par) -> ground_intensity(par[1], h, t)
    prob = QuadratureProblem(f, min_time(h), t, par)
    sol = solve(prob, HCubatureJL())
    return sol[1]
end

@doc raw"""
    logpdf(pp, h)

Compute the log probability density function for a temporal point process `pp` applied to history `h`:
```math
    \log f(h) = \sum_i \log \lambda(t_i | h_i) - \Lambda(h)
```

The default method uses a loop over events combined with [`integrated_ground_intensity`](@ref), but it should be reimplemented for specific processes if faster computation is possible.
"""
function Distributions.logpdf(pp::TemporalPointProcess, h::TemporalHistory)
    l = -integrated_ground_intensity(pp, h)
    for (t, m) in zip(event_times(h), event_marks(h))
        l += log(intensity(pp, h, t, m))
    end
    return l
end

@doc raw"""
    fit(pp0, h)

Compute the optimal parameter for a temporal point process of type `typeof(pp0)` on history `h` using maximum likelihood:
```math
    \hat{\theta} = \mathrm{argmax} \{f_{\theta}(h): \theta \in \Theta\}
```

The default method uses [GalacticOptim.jl](https://github.com/SciML/GalacticOptim.jl) for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.
"""
function Distributions.fit(
    pp_init::PP,
    h::TemporalHistory{M};
    adtype = GalacticOptim.AutoForwardDiff(),
    alg = Optim.LBFGS(),
) where {M,PP<:TemporalPointProcess}
    trans = build_transform(pp_init)
    θ_init = inverse(trans, ntfromstruct(pp_init))
    par = [nothing]
    f = (θ, par) -> -logpdf(PP(transform(trans, θ)), h)
    obj = OptimizationFunction(f, adtype)
    prob = OptimizationProblem(obj, θ_init, par)
    sol = solve(prob, alg)
    θ_opt = sol.u
    pp_opt = PP(transform(trans, θ_opt))
    return pp_opt
end

@doc raw"""
    check_residuals(pp, h)

Check whether the point process `pp` is a good fit for history `h` by applying Ogata's time rescaling method: if $(t_i)_i$ is a temporal point process with intensity $\lambda(t)$, then $(\Lambda(t_i))_i$ is a standard temporal Poisson process.
"""
function check_residuals(pp::TemporalPointProcess, h::TemporalHistory)
    Λ(t) = integrated_ground_intensity(pp, h, t)
    h_rescaled = time_change(h, Λ)
    qqplot_interevent_times(h_rescaled)
end
