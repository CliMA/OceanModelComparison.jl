# Instability of a perturbed Bickley jet

## Crude computational cost estimates:

For the "cost per DOF per time-step" estimate:

```
cost = simulation_time / (time_steps * DOF)
```

DOF      | Code                | Hardware                                      | Cost
:---:    | :---:               | :---:                                         | :---:
*32^2*   | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 4.81e-6
*32^2*   | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.64e-5
*32^2*   | Oceananigans.jl     | Titan V GPU                                   | 5.45e-5

*64^2*   | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 1.37e-6
*64^2*   | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 4.02e-6
*64^2*   | Oceananigans.jl     | Titan V GPU                                   | 7.94e-6

*128^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 4.67e-7
*128^2*  | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.49e-7
*128^2*  | Oceananigans.jl     | Titan V GPU                                   | 1.45e-6

*256^2*  | Oceananigans.jl     | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 1.20e-6
*256^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 7.93e-7
*256^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 6.98e-7
*256^2*  | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 5.15e-7
*256^2*  | Oceananigans.jl     | Titan V GPU                                   | 2.73e-7

*512^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 2.83e-7
*512^2*  | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 2.39e-7
*512^2*  | Oceananigans.jl     | Titan V GPU                                   | 5.83e-8

*1024^2* | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.91e-7
*1024^2* | Oceananigans.jl     | Titan V GPU                                   | 1.84e-8
