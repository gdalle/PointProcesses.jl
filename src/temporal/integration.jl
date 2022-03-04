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
    h::History,
    t = max_time(h),
)
    par = [pp]
    f = (t, par) -> ground_intensity(par[1], h, t)
    prob = QuadratureProblem(f, min_time(h), t, par)
    sol = solve(prob, HCubatureJL())
    return sol[1]
end
