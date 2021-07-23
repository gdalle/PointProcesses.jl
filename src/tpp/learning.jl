@doc raw"""
    integrated_ground_intensity(pp, h[, t])

Compute the integrated ground intensity (or compensator) $\Lambda(t|h)$ for a point process `pp` applied to history `h`:
```math
    \Lambda(h) = \int \lambda_g(t|h) \mathrm{d}t.
```

The default method uses the package [Quadrature](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit computation is feasible.
"""
function integrated_ground_intensity(
    pp::TemporalPointProcess{M},
    h::History{M},
    t::Float64 = h.tmax,
) where {M}
    prob = QuadratureProblem((t, pp) -> ground_intensity(pp, h, t), h.tmin, t, pp)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

@doc raw"""
    logpdf(pp, h)

Compute the log probability density function for a point process `pp` applied to history `h`:
```math
    \log f(h) = \sum_i \log \lambda(t_i | h_i) - \Lambda(h)
```

The default method uses the package [Quadrature](https://github.com/SciML/Quadrature.jl) for numerical integration, but it should be reimplemented for specific processes if explicit integration is feasible.
"""
function Distributions.logpdf(pp::TemporalPointProcess{M}, h::History{M}) where {M}
    l = -integrated_ground_intensity(pp, h)
    for (t, m) in zip(h.times, h.marks)
        l += log(intensity(pp, h, t, m))
    end
    return l
end

@doc raw"""
    fit(pp0, h)

Compute the optimal parameter for a point process of type `typeof(pp0)` on history `h` using maximum likelihood:
```math
    \hat{\theta} = \mathrm{argmax} \{f_{\theta}(h): \theta \in \Theta\}
```

The default method uses the package [GalacticOptim](https://github.com/SciML/GalacticOptim.jl) for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.
"""
function Distributions.fit(pp_init::PP, h::History{M}) where {M,PP<:TemporalPointProcess{M}}
    t = build_transform(pp_init)
    θ_init = inverse(t, params(pp_init))
    f = θ -> -logpdf(PP(transform(t, θ)), h)
    res = optimize(f, θ_init, LBFGS(), autodiff = :forward)
    θ_opt = Optim.minimizer(res)
    pp_opt = PP(transform(t, θ_opt))
    return pp_opt
end

@doc raw"""
    check_residuals(pp, h)

Check whether the point process `pp` is a good fit for history `h` by applying Ogata's time rescaling method: if $(t_i)_i$ is a point process with intensity $\lambda(t)$, then $(\Lambda(t_i))_i$ is a unit Poisson process.
"""
function check_residuals(pp::TemporalPointProcess{M}, h::History{M}) where {M}
    Λ(t) = integrated_ground_intensity(pp, h, t)
    h_rescaled = time_change(h, Λ)
    qqplot_interevent_times(h_rescaled)
end
