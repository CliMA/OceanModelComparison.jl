# Instability of a perturbed Bickley jet

## Crude computational cost estimates:

For the "cost per DOF per time-step" estimate:

```
cost = simulation_time / (time_steps * DOF)
```

DOF     | Code                | Hardware                                      | Cost
:---:   | :---:               | :---:                                         | :---:
*32^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 4.81e-6
*64^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 1.37e-6
*128^2* | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 4.67e-7
*256^2* | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 7.93e-7
*512^2* | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 7.93e-7


### Oceananigans.jl

### Exasim



We find


