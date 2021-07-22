@doc raw"""
    integrated_ground_intensity(pptype, θ, h)

Compute the integrated ground intensity (or compensator) $\Lambda(t|h)$ for a point process of type `pptype` with parameter `θ` applied to history `h`:
```math
    \Lambda(h) = \int \lambda_g(t|h) \mathrm{d}t.
```

The default method uses the package [Quadrature](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit computation is feasible.
"""
function integrated_ground_intensity(
    pptype::Type{<:PointProcess{M}},
    θ::Parameter,
    h::History{M},
) where {M}
    prob = QuadratureProblem((t, θ) -> ground_intensity(pptype, θ, h, t), h.tmin, h.tmax, θ)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

"""
    integrated_ground_intensity(pp, h)

Fall back on [`integrated_ground_intensity(pptype, θ, h)`](@ref).
"""
function integrated_ground_intensity(pp::PointProcess{M}, h::History{M}) where {M}
    return integrated_ground_intensity(typeof(pp), params(pp), h)
end

@doc raw"""
    logpdf(pptype, θ, h)

Compute the log probability density function for a point process of type `pptype` with parameter `θ` applied to history `h`:
```math
    \log f(h) = \sum_i \log \lambda(t_i | h_i) - \Lambda(h)
```

The default method uses the package [Quadrature](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.
"""
function Distributions.logpdf(
    pptype::Type{<:PointProcess{M}},
    θ::Parameter,
    h::History{M},
) where {M}
    l = -integrated_ground_intensity(pptype, θ, h)
    for (t, m) in zip(h.times, h.marks)
        l += log(intensity(pptype, θ, h, t, m))
    end
    return l
end

"""
    logpdf(pp, h)

Fall back on [`logpdf(pptype, θ, h)`](@ref).
"""
function Distributions.logpdf(pp::PointProcess{M}, h::History{M}) where {M}
    return logpdf(typeof(pp), params(pp), h)
end

@doc raw"""
    fit(pptype, h)

Compute the optimal parameter for a point process of type `pptype` on history `h` using maximum likelihood:
```math
    \hat{\theta} = \mathrm{argmax} \{f_{\theta}(h): \theta \in \Theta\}
```

The default method uses the package [GalacticOptim](https://github.com/SciML/GalacticOptim.jl) for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.
"""
function Distributions.fit(pptype::Type{<:PointProcess{M}}, θ0::Parameter, h::History{M}) where {M}
    f = OptimizationFunction(
        (θ, p) -> -logpdf(pptype, θ, h),
        GalacticOptim.AutoForwardDiff(),
    )
    prob = OptimizationProblem(f, θ0)
    sol = solve(prob, LBFGS())
    θ_opt = sol.minimizer
    pp_opt = pptype(NamedTuple(θ_opt)...)
    return pp_opt
end

function Distributions.fit(pp0::PointProcess{M}, h::History{M}) where {M}
    return fit(typeof(pp0), params(pp0), h)
end
