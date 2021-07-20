"""
    PointProcess{M}

The common supertype of all point processes with mark type `M`.
"""
abstract type PointProcess{M} end

"""
    Parameter

An alias for `ComponentVector`, which is how we store point process parameters `θ`.
"""
const Parameter = ComponentVector

"""
    get_θ(pp::PointProcess)

Build a `Parameter` from a `NamedTuple` linking the fields and values of `pp`.
"""
function get_θ(pp::PointProcess)
    return ComponentVector{Real}(ntfromstruct(pp))
end

# Functions with mandatory implementation by concrete types, and their fallbacks

@doc raw"""
    intensity(pptype, θ, h, t, m)

Compute the conditional intensity for a point process of type `pptype` with parameter `θ` applied to history `h` and event `(t, m)`.

The conditional intensity function quantifies the instantaneous risk  of an event with mark `m` occurring at time `t`[^Rasmussen_2018].

[^Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].
"""
function intensity(pptype, θ, h, t, m) end

function intensity(pp::PointProcess{M}, h::History{M}, t::Real, m::M) where {M}
    return intensity(typeof(pp), get_θ(pp), h, t, m)
end

@doc raw"""
    mark_distribution(pptype, θ, h, t)

Compute the distribution of marks for a point process of type `pptype` with parameter `θ` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(pptype, θ, h, t) end

function mark_distribution(pp::PointProcess{M}, h::History{M}, t::Real) where {M}
    return mark_distribution(typeof(pp), get_θ(pp), h, t)
end

@doc raw"""
    ground_intensity(pptype, θ, h, t)

Compute the ground intensity for a point process of type `pptype` with parameter `θ` applied to history `h` at time `t`.

The ground intensity quantifies the instantaneous risk of an event with any mark occurring at time `t`[^Rasmussen_2018]. It can be expressed as
```math
    \lambda_g(t|h) = \sum_{m \in \mathcal{M}} \lambda(t, m|h).
```
"""
function ground_intensity(pptype, θ, h, t) end

function ground_intensity(pp::PointProcess{M}, h::History{M}, t::Real) where {M}
    return ground_intensity(typeof(pp), get_θ(pp), h, t)
end

@doc raw"""
    ground_intensity_bound(pptype, θ, h, t)

Compute a local upper bound on the ground intensity for a point process of type `pptype` with parameter `θ` applied to history `h` at time `t`[^Rasmussen_2018].

Return a tuple of the form $(B, L)$ satisfying
```math
    \forall u \in [t, t+L), \quad \lambda_g(t|h) \leq B
```
"""
function ground_intensity_bound(pptype, θ, h, t) end

function ground_intensity_bound(pp::PointProcess{M}, h::History{M}, t::Real) where {M}
    return ground_intensity_bound(typeof(pp), get_θ(pp), h, t)
end

# Functions with optional implementation by concrete types

@doc raw"""
    integrated_ground_intensity(pptype, θ, h)

Compute the integrated ground intensity (or compensator) $\Lambda(t|h)$ for a point process of type `pptype` with parameter `θ` applied to history `h`:
```math
    \Lambda(h) = \int \lambda_g(t|h) \mathrm{d}t.
```

The default method uses the package [Quadrature](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit computation is feasible.
"""
function integrated_ground_intensity(pptype::Type{<:PointProcess}, θ, h)
    prob = QuadratureProblem((t, θ) -> ground_intensity(pptype, θ, h, t), h.tmin, h.tmax, θ)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

function integrated_ground_intensity(pp::PointProcess{M}, h::History{M}) where {M}
    return integrated_ground_intensity(typeof(pp), get_θ(pp), h)
end

@doc raw"""
    logpdf(pptype, θ, h)

Compute the log probability density function for a point process of type `pptype` with parameter `θ` applied to history `h`:
```math
    \log f(h) = \sum_i \log \lambda(t_i | h_i) - \Lambda(h)
```

The default method uses the package [Quadrature](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.
"""
function Distributions.logpdf(pptype::Type{<:PointProcess}, θ, h)
    l = -integrated_ground_intensity(pptype, θ, h)
    for (t, m) in zip(h.times, h.marks)
        l += log(intensity(pptype, θ, h, t, m))
    end
    return l
end

function Distributions.logpdf(pp::PointProcess{M}, h::History{M}) where {M}
    return logpdf(typeof(pp), get_θ(pp), h)
end

## Learning

@doc raw"""
    fit(pptype, h)

Compute the optimal parameter for a point process of type `pptype` on history `h` using maximum likelihood:
```math
    \hat{\theta} = \mathrm{argmax} \{f_{\theta}(h): \theta \in \Theta\}
```

The default method uses the package [GalacticOptim](https://github.com/SciML/GalacticOptim.jl) for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.

# Examples

```jldoctest
julia> using Random; Random.seed!(96);

julia> pp = MultivariatePoissonProcess([0., 1., 2.]);

julia> h = rand(pp, 0., 1000.);

julia> θ_est = fit(MultivariatePoissonProcess, h);

julia> pp_est = MultivariatePoissonProcess(θ_est.logλ)
MultivariatePoissonProcess([-0.008032171697071258, 0.9895411936136185, 2.0016151262152886])
```
"""
function Distributions.fit(pptype::Type{<:PointProcess}, h)
    f = OptimizationFunction(
        (θ, params) -> -logpdf(pptype, θ, h),
        GalacticOptim.AutoForwardDiff(),
    )
    prob = OptimizationProblem(f, default_θ(pptype, h))
    sol = solve(prob, LBFGS())
    θ_opt = sol.minimizer
    return θ_opt
end

# Simulation

function Base.rand(rng::AbstractRNG, pp::PointProcess{M}, tmin::Real, tmax::Real) where {M}
    h = History(Float64[], M[], tmin, tmax)
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, h, t)
        T = B > 0 ? rand(Exponential(1 / B)) : Inf
        if T > L
            t = t + L
        elseif T <= L
            U = rand(Uniform(0, 1))
            if U < ground_intensity(pp, h, t + T) / B
                m = rand(mark_distribution(pp, h, t + T))
                if t + T < tmax
                    push!(h, t + T, m)
                end
            end
            t = t + T
        end
    end
    return h
end

"""
    rand(pp, tmin, tmax)

Simulate a point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm[^Rasmussen_2018].

# Examples

```jldoctest
julia> pp = MultivariatePoissonProcess([0., 1., 2.]);

julia> h = rand(pp, 0., 1000.);
```
"""
Base.rand(pp::PointProcess, tmin::Real, tmax::Real) = rand(GLOBAL_RNG, pp, tmin, tmax)
