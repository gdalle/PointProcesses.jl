using PointProcesses
using Distributions

poissonprocess = PoissonProcess(1., Categorical([0.1, 0.9]))

history = simulate(poissonprocess; tmin=0, tmax=100)

plot(history)