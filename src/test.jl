using ParameterHandling # load up this package.
using Optim # generic optimisation
using Zygote # algorithmic differentiation
using AbstractGPs # package containing the models we'll be working with
using GalacticOptim

# Declare a NamedTuple containing an initial guess at parameters.
raw_initial_params = (
    k1 = (
        var=positive(0.9),
        precision=positive(1.0),
    ),
    k2 = (
        var=positive(0.1),
        precision=positive(0.3),
    ),
    noise_var=positive(0.2),
)

# Using ParameterHandling.flatten, we can obtain both a Vector{Float64} representation of
# these parameters, and a mapping from that vector back to the original parameters:
flat_initial_params, unflatten = ParameterHandling.flatten(raw_initial_params)

# ParameterHandling.value strips out all of the Positive types in initial_params,
# returning a plain named tuple of named tuples and Float64s.
initial_params = ParameterHandling.value(raw_initial_params)

# We define unpack to map directly from the flat Vector{Float64} representation to a
# the NamedTuple representation with all the Positive types removed.
unpack = ParameterHandling.value ∘ unflatten

# GP-specific functionality. Don't worry about the details, just
# note the use of the structured representation of the parameters.
function build_gp(params::NamedTuple)
    k1 = params.k1.var * Matern52Kernel() ∘ ScaleTransform(params.k1.precision)
    k2 = params.k2.var * SEKernel() ∘ ScaleTransform(params.k2.precision)
    return GP(k1 + k2)
end

# Generate some synthetic training data.
# Again, don't worry too much about the specifics here.
const x = range(-5.0, 5.0; length=100)
const y = rand(build_gp(initial_params)(x, initial_params.noise_var))

# Specify an objective function in terms of x and y.
function objective(params::NamedTuple)
    f = build_gp(params)
    return -logpdf(f(x, params.noise_var), y)
end

# Use Optim.jl to minimise the objective function w.r.t. the params.
# The important thing here is to note that we're passing in the flat vector of parameters to
# Optim, which is something that Optim knows how to work with, and using `unpack` to convert
# from this representation to the structured one that our objective function knows about
# using `unpack` -- we've used ParameterHandling to build a bridge between Optim and an
# entirely unrelated package.
training_results = Optim.optimize(
    objective ∘ unpack,
    θ -> only(Zygote.gradient(objective ∘ unpack, θ)),
    flat_initial_params,
    BFGS(
        alphaguess = Optim.LineSearches.InitialStatic(scaled=true),
        linesearch = Optim.LineSearches.BackTracking(),
    ),
    Optim.Options(show_trace = true);
    inplace=false,
)


f = OptimizationFunction((θ, p) -> (objective ∘ unpack)(θ), GalacticOptim.AutoZygote())
p = [1.]
prob = OptimizationProblem(f, flat_initial_params, p)
sol = solve(prob, BFGS())
x_opt = sol.minimizer

# Extracting the final values of the parameters.
final_params = unpack(training_results.minimizer)
f_trained = build_gp(final_params)
