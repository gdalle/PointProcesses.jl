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
    prob = QuadratureProblem((t, pp) -> ground_intensity(pp, h, t), h.tmin, t, pp)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
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

"""
    params(pp)

Retrieve point process parameters as a `NamedTuple`.
"""
StatsBase.params(pp::AbstractPointProcess) = ntfromstruct(pp)

@doc raw"""
    fit(pp0, h)

Compute the optimal parameter for a temporal point process of type `typeof(pp0)` on history `h` using maximum likelihood:
```math
    \hat{\theta} = \mathrm{argmax} \{f_{\theta}(h): \theta \in \Theta\}
```

The default method uses [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl/) for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.
"""
function Distributions.fit(pp_init::P, h::TemporalHistory{M}) where {M, P<:TemporalPointProcess{M}}
    t = build_transform(pp_init)
    θ_init = inverse(t, params(pp_init))
    f = θ -> -logpdf(P(transform(t, θ)), h)
    res = optimize(f, θ_init, LBFGS(), autodiff = :forward)
    θ_opt = Optim.minimizer(res)
    pp_opt = P(transform(t, θ_opt))
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
