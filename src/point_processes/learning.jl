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
    t = h.tmax,
)
    par = [pp]
    f = (t, par) -> ground_intensity(par[1], h, t)
    prob = QuadratureProblem(f, h.tmin, t, par)
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
    for (t, m) in zip(h.times, h.marks)
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
