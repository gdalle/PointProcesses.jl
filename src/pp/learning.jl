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

# Examples

```jldoctest
Random.seed!(96);
pp = MultivariatePoissonProcess([0., 1., 2.]);
h = rand(pp, 0., 1000.);
θ_est = fit(MultivariatePoissonProcess, h);
pp_est = MultivariatePoissonProcess(θ_est.logλ)

# output

MultivariatePoissonProcess([-0.008032171697071258, 0.9895411936136185, 2.0016151262152886])
```
"""
function Distributions.fit(pptype::Type{<:PointProcess{M}}, h::History{M}) where {M}
    f = OptimizationFunction(
        (θ, p) -> -logpdf(pptype, θ, h),
        GalacticOptim.AutoForwardDiff(),
    )
    prob = OptimizationProblem(f, default_param(pptype, h))
    sol = solve(prob, LBFGS())
    θ_opt = sol.minimizer
    return θ_opt
end
