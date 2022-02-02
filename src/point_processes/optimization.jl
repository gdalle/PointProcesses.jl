@doc raw"""
    fit(pp0, h)

Compute the optimal parameter for a temporal point process of type `typeof(pp0)` on history `h` using maximum likelihood:
```math
    \hat{\theta} = \mathrm{argmax} \{f_{\theta}(h): \theta \in \Theta\}
```

The default method uses [GalacticOptim.jl](https://github.com/SciML/GalacticOptim.jl) for numerical optimization, but it should be reimplemented for specific processes if explicit maximization is feasible.
"""
function fit(
    pp_init::PP,
    h::History{M};
    adtype = SciMLBase.NoAD(),
    alg = Optim.NelderMead(),
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
