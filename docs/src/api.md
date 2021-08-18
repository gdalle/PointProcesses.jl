# API reference

## Event history

### Abstract type

```@autodocs
Modules = [PointProcesses]
Pages = ["history/abstract.jl"]
```

### Temporal history

```@autodocs
Modules = [PointProcesses]
Pages = ["history/temporal.jl"]
```

## Point processes

### Abstract type

```@autodocs
Modules = [PointProcesses]
Pages = ["point_processes/abstract.jl"]
```

### Temporal point processes

```@autodocs
Modules = [PointProcesses]
Pages = ["point_processes/temporal.jl"]
```


## Built-in models

### Poisson processes

```@autodocs
Modules = [PointProcesses]
Pages = ["models/poisson.jl", "models/poisson_inhomogeneous.jl"]
```

```@autodocs
Modules = [PointProcesses]
Pages = ["models/poisson_multivariate.jl", "models/poisson_multivariate_naive.jl"]
```

### Hawkes processes

```@autodocs
Modules = [PointProcesses]
Pages = ["models/hawkes.jl"]
```

## Markov chains

### Abstract type

```@autodocs
Modules = [PointProcesses]
Pages = ["markov/abstract.jl"]
```

### Discrete time

```@autodocs
Modules = [PointProcesses]
Pages = ["markov/discrete_time.jl"]
```

### Continuous time

```@autodocs
Modules = [PointProcesses]
Pages = ["markov/continuous_time.jl"]
```

## Hidden Markov models

### Discrete time

```@autodocs
Modules = [PointProcesses]
Pages = ["hmm/hmm.jl", "hmm/baum_welch.jl"]
```

### Continuous time

```@autodocs
Modules = [PointProcesses]
Pages = ["hmm/mmpp.jl", "hmm/ryden.jl"]
```

## Utilities

```@autodocs
Modules = [PointProcesses]
Order = [:type, :function]
Pages = ["utils/plot.jl", "utils/utils.jl", "utils/categorical.jl", "utils/mvcategorical.jl", "utils.randvals.jl"]
```