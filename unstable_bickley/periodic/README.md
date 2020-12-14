# Instability of a perturbed Bickley jet

## Crude computational cost estimates:

For the "cost per DOF per time-step" estimate:

```
cost = simulation_time / (time_steps * DOF)
```

Note, `GeophysicalFlows.jl` used an 4-stage Runge-Kutta time-stepping method,
whereas `Oceananigans.jl` used an 3-stage Runge-Kutta time-stepping method.
A more fair comparison might multiply the `GeophysicalFlows.jl` results by 3/4.

DOF      | Code                | Hardware                                      | Cost    | Advection scheme |
:---:    | :---:               | :---:                                         | :---:   | :---: |
*32^2*   | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 4.81e-6 | |
*32^2*   | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 2.00e-5 | `WENO5` advection |
*32^2*   | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.11e-5 | `UpwindBiasedFifthOrder` advection |
*32^2*   | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.64e-5 | |
*32^2*   | GeophysicalFlows.jl | Titan V GPU                                   | 9.18e-7 | |
*32^2*   | Oceananigans.jl     | Titan V GPU                                   | 8.74e-6 | |

*64^2*   | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 1.37e-6 | |
*64^2*   | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 5.02e-6 | `WENO5` advection |
*64^2*   | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 4.07e-6 | `UpwindBiasedFifthOrder` advection |
*64^2*   | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 4.02e-6 | |
*64^2*   | GeophysicalFlows.jl | Titan V GPU                                   | 2.29e-7 | |
*64^2*   | Oceananigans.jl     | Titan V GPU                                   | 6.27e-6 | |

*128^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 4.67e-7 | |
*128^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.60e-6 | `WENO5` advection |
*128^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.37e-7 | `UpwindBiasedFifthOrder` advection |
*128^2*  | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.49e-7 | |
*128^2*  | GeophysicalFlows.jl | Titan V GPU                                   | 5.72e-8 | |
*128^2*  | Oceananigans.jl     | Titan V GPU                                   | 1.19e-6 | |

*256^2*  | Oceananigans.jl     | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 1.20e-6 | |
*256^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 7.93e-7 | |
*256^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 5.96e-7 | `WENO5` advection |
*256^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 4.97e-7 | `UpwindBiasedFifthOrder` advection |
*256^2*  | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 5.15e-7 | |
*256^2*  | Oceananigans.jl     | Titan V GPU                                   | 2.73e-7 | |
*256^2*  | GeophysicalFlows.jl | Titan V GPU                                   | 1.38e-8 | |

*512^2*  | GeophysicalFlows.jl | 8 threads with 3.1 GHz Intel i7 (MacBook Pro) | 2.83e-7 | |
*512^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 3.34e-7 | `WENO5` advection | 
*512^2*  | Oceananigans.jl     | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 4.34e-7 | `UpwindBiasedFifthOrder` advection | 
*512^2*  | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 2.39e-7 | |
*512^2*  | Oceananigans.jl     | Titan V GPU                                   | 4.94e-8 | |
*512^2*  | GeophysicalFlows.jl | Titan V GPU                                   | 6.83e-9 | |

*1024^2* | GeophysicalFlows.jl | 24 threads with 2.2 GHz Intel Xeon (x86_64)   | 1.91e-7 | |
*1024^2* | Oceananigans.jl     | Titan V GPU                                   | 1.89e-8 | |
*1024^2* | GeophysicalFlows.jl | Titan V GPU                                   | 7.07e-9 | |

*2048^2* | Oceananigans.jl     | Titan V GPU                                   | 1.47e-8 | |
*2048^2* | GeophysicalFlows.jl | Titan V GPU                                   | 7.07e-9 | |
