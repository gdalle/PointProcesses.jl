"""
    PointProcess{M}

The common supertype of all point processes with mark type `M`.
"""
abstract type PointProcess{M} end

# Parameters

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
    return ComponentVector{Float64}(ntfromstruct(pp))
end

# Functions with mandatory implementation by concrete types, and their fallbacks

@doc raw"""
    intensity(pptype, θ, h, t, m)

Compute the conditional intensity for a point process of type `pptype` with parameter `θ` applied to history `h` and event `(t, m)`.

The conditional intensity function quantifies the instantaneous risk  of an event with mark `m` occurring at time `t`[^Rasmussen_2018].

[^Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].
"""
function intensity(
    pptype::Type{<:PointProcess{M}},
    θ::Parameter,
    h::History{M},
    t::Float64,
    m::M,
) where {M} end

function intensity(pp::PointProcess{M}, h::History{M}, t::Float64, m::M) where {M}
    return intensity(typeof(pp), get_θ(pp), h, t, m)
end

@doc raw"""
    mark_distribution(pptype, θ, h, t)

Compute the distribution of marks for a point process of type `pptype` with parameter `θ` knowing that an event takes place at time `t` after history `h`.
"""
function mark_distribution(
    pptype::Type{<:PointProcess{M}},
    θ::Parameter,
    h::History{M},
    t::Float64,
) where {M} end

function mark_distribution(pp::PointProcess{M}, h::History{M}, t) where {M}
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
function ground_intensity(
    pptype::Type{<:PointProcess{M}},
    θ::Parameter,
    h::History{M},
    t::Float64,
) where {M} end

function ground_intensity(pp::PointProcess{M}, h::History{M}, t) where {M}
    return ground_intensity(typeof(pp), get_θ(pp), h, t)
end

@doc raw"""
    ground_intensity_bound(pptype, θ, h, t)

Compute a local upper bound on the ground intensity for a point process of type `pptype` with parameter `θ` applied to history `h` at time `t`[^Daley_2003].

Return a tuple of the form $(B, L)$ satisfying
```math
    \forall u \in [t, t+L), \quad \lambda_g(t|h) \leq B
```

[^Daley_2003]: Daley, D. J., and Vere-Jones, D. (eds.) (2003), “7. Conditional Intensities and Likelihoods,” in An Introduction to the Theory of Point Processes: Volume I: Elementary Theory and Methods, Probability and its Applications, New York, NY: Springer, pp. 211–287. https://doi.org/10.1007/0-387-21564-6_7.
"""
function ground_intensity_bound(
    pptype::Type{<:PointProcess{M}},
    θ::Parameter,
    h::History{M},
    t::Float64,
) where {M} end

function ground_intensity_bound(pp::PointProcess{M}, h::History{M}, t) where {M}
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
