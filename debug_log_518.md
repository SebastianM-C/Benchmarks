Using
```julia
using StaticArrays
using Profile
using BenchmarkTools
using OrdinaryDiffEq

function lorenz(u,p,t)
  du1 = 10.0*(u[2]-u[1])
  du2 = u[1]*(28.0-u[3]) - u[2]
  du3 = u[1]*u[2] - (8/3)*u[3]

  return SVector{3}(du1,du2,du3)
end

function lorenz2(u,p,t)
    du1 = 10.0*(u[2]-u[1])
    du2 = u[1]*(28.0-u[3]) - u[2]
    du3 = u[1]*u[2] - (8/3)*u[3]

    du4 = 10.0*(u[2+3]-u[1+3])
    du5 = u[1+3]*(28.0-u[3+3]) - u[2+3]
    du6 = u[1+3]*u[2+3] - (8/3)*u[3+3]

    return SVector{6}(du1,du2,du3,du4,du5,du6)
end

const u0 = @SVector [1.0,0.0,0.0]
const u = vcat(u0,u0)
const tspan = (0.0,10.0)
const prob = ODEProblem(lorenz,u0,tspan)
# const prob = ODEProblem(lorenz2,u,tspan)

@time solve(prob, Tsit5(), save_everystep=false);

@time solve(prob, Tsit5(), save_everystep=false);
Profile.clear_malloc_data()
@timev solve(prob, Tsit5(), save_everystep=false);
exit()
```

I ran the script with `julia --track-allocation=user --inline=no test.jl` with the result
```
 & "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user --inline=no .\test.jl
  7.693043 seconds (6.12 M allocations: 336.579 MiB, 2.85% gc time)
  0.594925 seconds (934.62 k allocations: 45.458 MiB, 0.78% gc time)
  0.601419 seconds (934.62 k allocations: 45.458 MiB, 0.48% gc time)
elapsed time (ns): 601419499
gc time (ns):      2899100
bytes allocated:   47666264
pool allocs:       934620
GC pauses:         1
```

The `.mem` file for the script indicates allocations in the user function
```
        - function lorenz2(u,p,t)
    39744     du1 = 10.0*(u[2]-u[1])
    59616     du2 = u[1]*(28.0-u[3]) - u[2]
    69552     du3 = u[1]*u[2] - (8/3)*u[3]
        -
    39744     du4 = 10.0*(u[2+3]-u[1+3])
    59616     du5 = u[1+3]*(28.0-u[3+3]) - u[2+3]
    69552     du6 = u[1+3]*u[2+3] - (8/3)*u[3+3]
        -
        0     return SVector{6}(du1,du2,du3,du4,du5,du6)
        - end
```

Looking in the directory of `OrdinaryDiffEq` with
```julia
using Coverage, OrdinaryDiffEq
info = analyze_malloc(dirname(pathof(OrdinaryDiffEq)))
info = info[findall(x->x.bytes != 0, info)]
```
shows allocations in `perform_step`.

First, in solve:
```
        - function DiffEqBase.solve!(integrator::ODEIntegrator)
        0   @inbounds while !isempty(integrator.opts.tstops)
     6656     while integrator.tdir * integrator.t < top(integrator.opts.tstops)
        0       loopheader!(integrator)
        0       if check_error!(integrator) != :Success
        0         return integrator.sol
        -       end
    57680       perform_step!(integrator,integrator.cache)
        0       loopfooter!(integrator)
        0       if isempty(integrator.opts.tstops)
        -         break
        -       end
        -     end
        0     handle_tstop!(integrator)
        -   end
        0   postamble!(integrator)
        -
       96   f = integrator.sol.prob.f
        -
        0   if DiffEqBase.has_analytic(f)
        -     DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
        -   end
        0   if integrator.sol.retcode != :Default
        0     return integrator.sol
        -   end
        0   integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
        - end
```
and then in `perform_step/low_order_rk_perform_step.jl`
```
        - @muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step=false)
    75808   @unpack t,dt,uprev,u,f,p = integrator
    56032   @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
     6592   k1 = integrator.fsalfirst
     1648   a = dt*a21
    14832   k2 = f(uprev+a*k1, p, t+c1*dt)
    34608   k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
    47792   k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
    60976   k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
    72512   g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
     1648   k6 = f(g6, p, t+dt)
    85696   u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    14832   integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
     1200   integrator.destats.nf += 6
        0   if typeof(integrator.alg) <: CompositeAlgorithm
        -     g7 = u
        -     # Hairer II, page 22
        -     integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
        -   end
        0   if integrator.opts.adaptive
    92288     utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
     3296     atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
     4944     integrator.EEst = integrator.opts.internalnorm(atmp,t)
        -   end
     6592   integrator.k[1] = k1
     6592   integrator.k[2] = k2
     6592   integrator.k[3] = k3
     6592   integrator.k[4] = k4
     6592   integrator.k[5] = k5
     6592   integrator.k[6] = k6
     6592   integrator.k[7] = k7
     6592   integrator.u = u
        - end
```

Now to compare with the case that doesn't present the problem I used
```
const prob = ODEProblem(lorenz,u0,tspan)
# const prob = ODEProblem(lorenz2,u,tspan)
```
and deleted all the `*.mem` files in `.julia` and ran again.
```
& "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user --inline=no .\test.jl
  8.038784 seconds (5.68 M allocations: 304.809 MiB, 8.82% gc time)
  0.564761 seconds (667.71 k allocations: 22.686 MiB, 1.68% gc time)
  0.594728 seconds (667.71 k allocations: 22.686 MiB)
elapsed time (ns): 594727700
bytes allocated:   23787592
pool allocs:       667713
```

Now the user function has
```
        - function lorenz(u,p,t)
    39744   du1 = 10.0*(u[2]-u[1])
    59616   du2 = u[1]*(28.0-u[3]) - u[2]
    69552   du3 = u[1]*u[2] - (8/3)*u[3]
        -
        0   return SVector{3}(du1,du2,du3)
        - end
```
solve
```
function DiffEqBase.solve!(integrator::ODEIntegrator)
        0   @inbounds while !isempty(integrator.opts.tstops)
     6656     while integrator.tdir * integrator.t < top(integrator.opts.tstops)
        0       loopheader!(integrator)
        0       if check_error!(integrator) != :Success
        0         return integrator.sol
        -       end
    54384       perform_step!(integrator,integrator.cache)
        0       loopfooter!(integrator)
        0       if isempty(integrator.opts.tstops)
        -         break
        -       end
        -     end
        0     handle_tstop!(integrator)
        -   end
        0   postamble!(integrator)
        -
       80   f = integrator.sol.prob.f
        -
        0   if DiffEqBase.has_analytic(f)
        -     DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
        -   end
        0   if integrator.sol.retcode != :Default
        0     return integrator.sol
        -   end
        0   integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
        - end
 ```
 and perform_setep
 ```
         - @muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step=false)
    65920   @unpack t,dt,uprev,u,f,p = integrator
    56032   @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
     3296   k1 = integrator.fsalfirst
     1648   a = dt*a21
     8240   k2 = f(uprev+a*k1, p, t+c1*dt)
    18128   k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
    24720   k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
    31312   k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
    36256   g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
     1648   k6 = f(g6, p, t+dt)
    42848   u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
     8240   integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
     1200   integrator.destats.nf += 6
        0   if typeof(integrator.alg) <: CompositeAlgorithm
        -     g7 = u
        -     # Hairer II, page 22
        -     integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
        -   end
        0   if integrator.opts.adaptive
    46144     utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
     3296     atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
     4944     integrator.EEst = integrator.opts.internalnorm(atmp,t)
        -   end
     3296   integrator.k[1] = k1
     3296   integrator.k[2] = k2
     3296   integrator.k[3] = k3
     3296   integrator.k[4] = k4
     3296   integrator.k[5] = k5
     3296   integrator.k[6] = k6
     3296   integrator.k[7] = k7
     3296   integrator.u = u
        - end
```

Now, as a sanity check we double the time and run again for `lorenz` to see that the problem is not present
```
& "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user --inline=no .\test.jl
  8.996536 seconds (6.48 M allocations: 331.773 MiB, 7.13% gc time)
  1.259313 seconds (1.46 M allocations: 49.645 MiB, 0.52% gc time)
  1.238782 seconds (1.46 M allocations: 49.645 MiB, 0.37% gc time)
elapsed time (ns): 1238781701
gc time (ns):      4545901
bytes allocated:   52056944
pool allocs:       1460701
GC pauses:         1
```
This is quite strange, allocations increased. Let's try without disabling inlining.
For
```julia
const tspan = (0.0,10.0)
const prob = ODEProblem(lorenz,u0,tspan)
# const prob = ODEProblem(lorenz2,u,tspan)
```
we get
```
 & "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user .\test.jl
  8.964144 seconds (7.13 M allocations: 366.101 MiB, 3.86% gc time)
  0.002369 seconds (50 allocations: 4.719 KiB)
  0.002053 seconds (50 allocations: 4.719 KiB)
elapsed time (ns): 2053300
bytes allocated:   4832
pool allocs:       50
```
Now we don't have allocations during `perform_step~` in `solve!`
```
        - function DiffEqBase.solve!(integrator::ODEIntegrator)
        0   @inbounds while !isempty(integrator.opts.tstops)
        0     while integrator.tdir * integrator.t < top(integrator.opts.tstops)
        0       loopheader!(integrator)
        0       if check_error!(integrator) != :Success
        0         return integrator.sol
        -       end
        0       perform_step!(integrator,integrator.cache)
        0       loopfooter!(integrator)
        0       if isempty(integrator.opts.tstops)
        -         break
        -       end
        -     end
        0     handle_tstop!(integrator)
        -   end
        0   postamble!(integrator)
        -
        -   f = integrator.sol.prob.f
        -
        -   if DiffEqBase.has_analytic(f)
        -     DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
        -   end
        0   if integrator.sol.retcode != :Default
        0     return integrator.sol
        -   end
      128   integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
        - end
```
and increasing the time
```julia
const tspan = (0.0,20.0)
```
gives
```
& "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user .\test.jl
  7.858750 seconds (7.13 M allocations: 366.101 MiB, 4.89% gc time)
  0.004894 seconds (50 allocations: 4.719 KiB)
  0.005789 seconds (50 allocations: 4.719 KiB)
elapsed time (ns): 5789000
bytes allocated:   4832
pool allocs:       50
```
which has the same number of allocations as in the previous case. Thus does not present the problem.

Let's now try again for `lorenz2`.
```julia
const tspan = (0.0,10.0)
# const prob = ODEProblem(lorenz,u0,tspan)
const prob = ODEProblem(lorenz2,u,tspan)
```
The first run gives
```
& "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user .\test.jl
  8.975788 seconds (7.52 M allocations: 381.617 MiB, 4.78% gc time)
  0.002314 seconds (368 allocations: 20.156 KiB)
  0.003371 seconds (368 allocations: 20.156 KiB)
elapsed time (ns): 3371000
bytes allocated:   20640
pool allocs:       368
```

Using the Coverage script
```julia
using Coverage, OrdinaryDiffEq
info = analyze_malloc(dirname(pathof(OrdinaryDiffEq)))
info = info[findall(x->x.bytes != 0, info)]
```
we obtain
```
julia> info = info[findall(x->x.bytes != 0, info)]
24-element Array{Coverage.MallocInfo,1}:
 Coverage.MallocInfo(16, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\initdt.jl.3152.mem", 131)
 Coverage.MallocInfo(16, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\initdt.jl.3152.mem", 138)
 Coverage.MallocInfo(16, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\initdt.jl.3152.mem", 156)
 Coverage.MallocInfo(32, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 463)  
 Coverage.MallocInfo(48, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 214)  
 Coverage.MallocInfo(48, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 215)  
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 116)  
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 195)  
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 196)  
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 240)  
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 287)  
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 420)  
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 438)  
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 457)  
 Coverage.MallocInfo(128, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 224)
 Coverage.MallocInfo(144, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 295)
 Coverage.MallocInfo(144, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 390)
 Coverage.MallocInfo(208, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 213)
 Coverage.MallocInfo(224, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 267)
 Coverage.MallocInfo(448, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\perform_step\\low_order_rk_perform_step.jl.3152.mem", 564)
 Coverage.MallocInfo(496, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 262)
 Coverage.MallocInfo(832, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 66)  
 Coverage.MallocInfo(1008, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.3152.mem", 339)
 Coverage.MallocInfo(1648, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\perform_step\\low_order_rk_perform_step.jl.3152.mem", 599)
```

The largest allocation is in `perform_step!` on line 599 for the error estimator
```
- @muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step=false)
0   @unpack t,dt,uprev,u,f,p = integrator
0   @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
0   k1 = integrator.fsalfirst
0   a = dt*a21
0   k2 = f(uprev+a*k1, p, t+c1*dt)
0   k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
0   k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
0   k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
0   g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
0   k6 = f(g6, p, t+dt)
0   u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
0   integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
0   integrator.destats.nf += 6
-   if typeof(integrator.alg) <: CompositeAlgorithm
-     g7 = u
-     # Hairer II, page 22
-     integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
-   end
0   if integrator.opts.adaptive
0     utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
0     atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
1648     integrator.EEst = integrator.opts.internalnorm(atmp,t)
-   end
0   integrator.k[1] = k1
0   integrator.k[2] = k2
0   integrator.k[3] = k3
0   integrator.k[4] = k4
0   integrator.k[5] = k5
0   integrator.k[6] = k6
0   integrator.k[7] = k7
0   integrator.u = u
- end
```
It is a bit strange that if we look in `solve!`, `perform_step!` doesn't show
any allocations. (Can this be due to inlining?)
```
        - function DiffEqBase.solve!(integrator::ODEIntegrator)
        0   @inbounds while !isempty(integrator.opts.tstops)
        0     while integrator.tdir * integrator.t < top(integrator.opts.tstops)
        0       loopheader!(integrator)
        0       if check_error!(integrator) != :Success
        0         return integrator.sol
        -       end
        0       perform_step!(integrator,integrator.cache)
        0       loopfooter!(integrator)
        0       if isempty(integrator.opts.tstops)
        -         break
        -       end
        -     end
        0     handle_tstop!(integrator)
        -   end
        0   postamble!(integrator)
        -
        -   f = integrator.sol.prob.f
        -
        -   if DiffEqBase.has_analytic(f)
        -     DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
        -   end
        0   if integrator.sol.retcode != :Default
        0     return integrator.sol
        -   end
      144   integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
        - end
```

The second largest allocation is in `__init` on line 339 for the integrator.
```
- function DiffEqBase.__init(prob::Union{DiffEqBase.AbstractODEProblem,DiffEqBase.AbstractDAEProblem},
-                            alg::Union{OrdinaryDiffEqAlgorithm,DAEAlgorithm},
-                            timeseries_init = typeof(prob.u0)[],
-                            ts_init = eltype(prob.tspan)[],
-                            ks_init = [],
-                            recompile::Type{Val{recompile_flag}} = Val{true};
-                            saveat = eltype(prob.tspan)[],
-                            tstops = eltype(prob.tspan)[],
-                            d_discontinuities= eltype(prob.tspan)[],
-                            save_idxs = nothing,
-                            save_everystep = isempty(saveat),
-                            save_on = true,
-                            save_start = save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[1] in saveat,
-                            save_end = save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[2] in saveat,
-                            callback = nothing,
-                            dense = save_everystep && !(typeof(alg) <: Union{DAEAlgorithm,FunctionMap}) && isempty(saveat),
-                            calck = (callback !== nothing && callback != CallbackSet()) || # Empty callback
-                                    (!isempty(setdiff(saveat,tstops)) || dense), # and no dense output
-                            dt = typeof(alg) <: FunctionMap && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
-                            dtmin = typeof(one(eltype(prob.tspan))) <: AbstractFloat ? eps(eltype(prob.tspan)) :
-                                    typeof(one(eltype(prob.tspan))) <: Integer ? 0 : eltype(prob.tspan)(1//10^(10)),
-                            dtmax = eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
-                            force_dtmin = false,
-                            adaptive = isadaptive(alg),
-                            gamma = gamma_default(alg),
-                            abstol = nothing,
-                            reltol = nothing,
-                            qmin = qmin_default(alg),
-                            qmax = qmax_default(alg),
-                            qsteady_min = qsteady_min_default(alg),
-                            qsteady_max = qsteady_max_default(alg),
-                            qoldinit = 1//10^4,
-                            fullnormalize = true,
-                            failfactor = 2,
-                            beta1 = nothing,
-                            beta2 = nothing,
-                            maxiters = adaptive ? 1000000 : typemax(Int),
-                            internalnorm = ODE_DEFAULT_NORM,
-                            internalopnorm = LinearAlgebra.opnorm,
-                            isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
-                            unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
-                            verbose = true,
-                            timeseries_errors = true,
-                            dense_errors = false,
-                            advance_to_tstop = false,
-                            stop_at_next_tstop = false,
-                            initialize_save = true,
-                            progress = false,
-                            progress_steps = 1000,
-                            progress_name = "ODE",
-                            progress_message = ODE_DEFAULT_PROG_MESSAGE,
-                            userdata = nothing,
-                            allow_extrapolation = alg_extrapolates(alg),
-                            initialize_integrator = true,
-                            alias_u0 = false,
-                            kwargs...) where recompile_flag
-
832   if prob isa DiffEqBase.AbstractDAEProblem && alg isa OrdinaryDiffEqAlgorithm
-     error("You cannot use an ODE Algorithm with a DAEProblem")
-   end
-
-   if prob isa DiffEqBase.AbstractODEProblem && alg isa DAEAlgorithm
-     error("You cannot use an DAE Algorithm with a ODEProblem")
-   end
-
-   if typeof(prob.f)<:DynamicalODEFunction && typeof(prob.f.mass_matrix)<:Tuple
-     if any(mm != I for mm in prob.f.mass_matrix)
-       error("This solver is not able to use mass matrices.")
-     end
0   elseif !(typeof(prob)<:DiscreteProblem) &&
-          !(typeof(prob)<:DiffEqBase.AbstractDAEProblem) &&
-          !is_mass_matrix_alg(alg) &&
-          prob.f.mass_matrix != I
0     error("This solver is not able to use mass matrices.")
-   end
-
0   if !isempty(saveat) && dense
0     @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
-   end
-
0   progress && @logmsg(LogLevel(-1),progress_name,_id=_id = :OrdinaryDiffEq,progress=0)
-
0   tType = eltype(prob.tspan)
0   tspan = prob.tspan
0   tdir = sign(tspan[end]-tspan[1])
-
0   t = tspan[1]
-
0   if (((!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm) && !(typeof(alg) <: DAEAlgorithm)) || !adaptive) && dt == tType(0) && isempty(tstops)) && !(typeof(alg) <: Union{FunctionMap,LinearExponential})
0       error("Fixed timestep methods require a choice of dt or choosing the tstops")
-   end
-
0   f = prob.f
-   p = prob.p
-
-   # Get the control variables
-
0   if alias_u0
0     u = prob.u0
-   else
0     u = recursivecopy(prob.u0)
-   end
-
0   uType = typeof(u)
-   uBottomEltype = recursive_bottom_eltype(u)
-   uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)
-
80   ks = Vector{uType}(undef, 0)
-
-   uEltypeNoUnits = recursive_unitless_eltype(u)
-   tTypeNoUnits   = typeof(one(tType))
-
-   if typeof(alg) <: FunctionMap
-     abstol_internal = real.(zero.(u))
-   elseif abstol === nothing
-     if uBottomEltypeNoUnits == uBottomEltype
0       abstol_internal = real(convert(uBottomEltype,oneunit(uBottomEltype)*1//10^6))
-     else
-       abstol_internal = real.(oneunit.(u).*1//10^6)
-     end
-   else
-     abstol_internal = real.(abstol)
-   end
-
-   if typeof(alg) <: FunctionMap
-     reltol_internal = real.(zero(first(u)/t))
-   elseif reltol === nothing
-     if uBottomEltypeNoUnits == uBottomEltype
0       reltol_internal = real(convert(uBottomEltype,oneunit(uBottomEltype)*1//10^3))
-     else
-       reltol_internal = real.(oneunit.(u).*1//10^3)
-     end
-   else
-     reltol_internal = real.(reltol)
-   end
-
0   dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
-   # dtmin is all abs => does not care about sign already.
-
0   if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number && uBottomEltypeNoUnits == uBottomEltype # Could this be more efficient for other arrays?
-     if !(typeof(u) <: ArrayPartition)
-       rate_prototype = recursivecopy(u)
-     else
-       rate_prototype = similar(u, typeof.(oneunit.(recursive_bottom_eltype.(u.x))./oneunit(tType))...)
-     end
-   else
-     if uBottomEltypeNoUnits == uBottomEltype
-       rate_prototype = u
-     else # has units!
-       rate_prototype = u/oneunit(tType)
-     end
-   end
-   rateType = typeof(rate_prototype) ## Can be different if united
-
0   tstops_internal, saveat_internal, d_discontinuities_internal =
-     tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan)
-
-   callbacks_internal = CallbackSet(callback)
-
0   max_len_cb = DiffEqBase.max_vector_callback_length(callbacks_internal)
-   if max_len_cb isa VectorContinuousCallback
-     if isinplace(prob)
-       callback_cache = DiffEqBase.CallbackCache(u,max_len_cb.len,uBottomEltype,uBottomEltype)
-     else
-       callback_cache = DiffEqBase.CallbackCache(max_len_cb.len,uBottomEltype,uBottomEltype)
-     end
-   else
-     callback_cache = nothing
-   end
-
-   ### Algorithm-specific defaults ###
-   if save_idxs === nothing
-     ksEltype = Vector{rateType}
-   else
-     ks_prototype = rate_prototype[save_idxs]
-     ksEltype = Vector{typeof(ks_prototype)}
-   end
-
-   # Have to convert incase passed in wrong.
-   if save_idxs === nothing
-     timeseries = convert(Vector{uType},timeseries_init)
-   else
-     u_initial = u[save_idxs]
-     timeseries = convert(Vector{typeof(u_initial)},timeseries_init)
-   end
-   ts = convert(Vector{tType},ts_init)
80   ks = convert(Vector{ksEltype},ks_init)
80   alg_choice = Int[]
-
0   if !adaptive && save_everystep && tspan[2]-tspan[1] != Inf
0     dt == 0 ? steps = length(tstops) :
-               steps = ceil(Int,internalnorm((tspan[2]-tspan[1])/dt,tspan[1]))
0     sizehint!(timeseries,steps+1)
0     sizehint!(ts,steps+1)
0     sizehint!(ks,steps+1)
0   elseif save_everystep
0     sizehint!(timeseries,50)
0     sizehint!(ts,50)
0     sizehint!(ks,50)
0   elseif !isempty(saveat_internal)
0     sizehint!(timeseries,length(saveat_internal)+1)
0     sizehint!(ts,length(saveat_internal)+1)
0     sizehint!(ks,length(saveat_internal)+1)
-   else
208     sizehint!(timeseries,2)
48     sizehint!(ts,2)
48     sizehint!(ks,2)
-   end
-
0   if save_start
-     saveiter = 1 # Starts at 1 so first save is at 2
-     saveiter_dense = 1
0     copyat_or_push!(ts,1,t)
-     if save_idxs === nothing
0       copyat_or_push!(timeseries,1,u)
128       copyat_or_push!(ks,1,[rate_prototype])
-     else
-       copyat_or_push!(timeseries,1,u_initial,Val{false})
-       copyat_or_push!(ks,1,[ks_prototype])
-     end
-
-     if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
-       copyat_or_push!(alg_choice,1,1)
-     end
-   else
-     saveiter = 0 # Starts at 0 so first save is at 1
-     saveiter_dense = 0
-   end
-
0   QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits
-
80   k = rateType[]
-
-   if uses_uprev(alg, adaptive) || calck
0     uprev = recursivecopy(u)
-   else
-     # Some algorithms do not use `uprev` explicitly. In that case, we can save
-     # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
-     uprev = u
-   end
0   if allow_extrapolation
0     uprev2 = recursivecopy(u)
-   else
-     uprev2 = uprev
-   end
-
0   cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol_internal,p,calck,Val(isinplace(prob)))
-
-   if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
-     id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,dense,cache)
-     beta2 === nothing && ( beta2=beta2_default(alg.algs[cache.current]) )
-     beta1 === nothing && ( beta1=beta1_default(alg.algs[cache.current],beta2) )
-   else
496     id = InterpolationData(f,timeseries,ts,ks,dense,cache)
0     beta2 === nothing && ( beta2=beta2_default(alg) )
0     beta1 === nothing && ( beta1=beta1_default(alg,beta2) )
-   end
-
224   opts = DEOptions{typeof(abstol_internal),typeof(reltol_internal),QT,tType,
-                    typeof(internalnorm),typeof(internalopnorm),typeof(callbacks_internal),typeof(isoutofdomain),
-                    typeof(progress_message),typeof(unstable_check),typeof(tstops_internal),
-                    typeof(d_discontinuities_internal),typeof(userdata),typeof(save_idxs),
-                    typeof(maxiters),typeof(tstops),typeof(saveat),
-                    typeof(d_discontinuities)}(
-                        maxiters,save_everystep,adaptive,abstol_internal,
-                        reltol_internal,QT(gamma),QT(qmax),
-                        QT(qmin),QT(qsteady_max),
-                        QT(qsteady_min),QT(failfactor),tType(dtmax),
-                        tType(dtmin),internalnorm,internalopnorm,save_idxs,tstops_internal,saveat_internal,
-                        d_discontinuities_internal,
-                        tstops,saveat,d_discontinuities,
-                        userdata,progress,progress_steps,
-                        progress_name,progress_message,timeseries_errors,dense_errors,
-                        QT(beta1),QT(beta2),QT(qoldinit),dense,
-                        save_on,save_start,save_end,callbacks_internal,isoutofdomain,
-                        unstable_check,verbose,
-                        calck,force_dtmin,advance_to_tstop,stop_at_next_tstop)
-
96   destats = DiffEqBase.DEStats(0)
-
-   if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
-     sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,
-                       dense=dense,k=ks,interp=id,
-                       alg_choice=alg_choice,
-                       calculate_error = false, destats=destats)
-   else
144     sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,
-                       dense=dense,k=ks,interp=id,
-                       calculate_error = false, destats=destats)
-   end
-
-   if recompile_flag == true
-     FType = typeof(f)
-     SolType = typeof(sol)
-     cacheType = typeof(cache)
-   else
-     FType = Function
-     if alg isa OrdinaryDiffEqAlgorithm
-       SolType = DiffEqBase.AbstractODESolution
-       cacheType =  OrdinaryDiffEqCache
-     else
-       SolType = DiffEqBase.AbstractDAESolution
-       cacheType =  DAECache
-     end
-   end
-
-   # rate/state = (state/time)/state = 1/t units, internalnorm drops units
-   eigen_est = one(uBottomEltypeNoUnits)/one(tType)
-   tprev = t
-   dtcache = tType(dt)
-   dtpropose = tType(dt)
-   iter = 0
-   kshortsize = 0
-   reeval_fsal = false
-   u_modified = false
-   EEst = tTypeNoUnits(1)
-   just_hit_tstop = false
-   isout = false
-   accept_step = false
-   force_stepfail = false
-   last_stepfail = false
-   event_last_time = 0
-   vector_event_last_time = 1
-   last_event_error = zero(uBottomEltypeNoUnits)
-   dtchangeable = isdtchangeable(alg)
-   q11 = tTypeNoUnits(1)
-   success_iter = 0
-   erracc = tTypeNoUnits(1)
-   dtacc = tType(1)
-
1008   integrator = ODEIntegrator{typeof(alg),isinplace(prob),uType,tType,typeof(p),typeof(eigen_est),
-                              QT,typeof(tdir),typeof(k),SolType,
-                              FType,cacheType,
-                              typeof(opts),fsal_typeof(alg,rate_prototype),
-                              typeof(last_event_error),typeof(callback_cache)}(
-                              sol,u,k,t,tType(dt),f,p,uprev,uprev2,tprev,
-                              alg,dtcache,dtchangeable,
-                              dtpropose,tdir,eigen_est,EEst,QT(qoldinit),q11,
-                              erracc,dtacc,success_iter,
-                              iter,saveiter,saveiter_dense,cache,callback_cache,
-                              kshortsize,force_stepfail,last_stepfail,
-                              just_hit_tstop,event_last_time,vector_event_last_time,last_event_error,
-                              accept_step,
-                              isout,reeval_fsal,
-                              u_modified,opts,destats)
0   if initialize_integrator
0     initialize_callbacks!(integrator, initialize_save)
0     initialize!(integrator,integrator.cache)
0     save_start && typeof(alg) <: CompositeAlgorithm && copyat_or_push!(alg_choice,1,integrator.cache.current)
-   end
-
0   handle_dt!(integrator)
-
0   integrator
- end
```

It is followed by line 66 (probably for the creation of `prob`) in the same
function. After that follows line 262 for the interpolation data (also in the
same function).

The next is on line 564 in the `initialize` function
```
- function initialize!(integrator, cache::Tsit5ConstantCache)
0   integrator.kshortsize = 7
448   integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
0   integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
0   integrator.destats.nf += 1
-
-   # Avoid undefined entries if k is an array of arrays
0   integrator.fsallast = zero(integrator.fsalfirst)
0   integrator.k[1] = integrator.fsalfirst
0   @inbounds for i in 2:integrator.kshortsize-1
0     integrator.k[i] = zero(integrator.fsalfirst)
-   end
0   integrator.k[integrator.kshortsize] = integrator.fsallast
- end
```


Let us now consider
```julia
const tspan = (0.0,20.0)
```
This gives
```
& "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user .\test.jl
  9.308793 seconds (7.52 M allocations: 381.634 MiB, 5.08% gc time)
  0.003696 seconds (737 allocations: 37.453 KiB)
  0.006386 seconds (737 allocations: 37.453 KiB)
elapsed time (ns): 6385600
bytes allocated:   38352
pool allocs:       737
```
which has roughly double the allocations. Analyzing with the Coverage script gives
```julia
julia> info = info[findall(x->x.bytes != 0, info)]
24-element Array{Coverage.MallocInfo,1}:
 Coverage.MallocInfo(16, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\initdt.jl.16956.mem", 131)
 Coverage.MallocInfo(16, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\initdt.jl.16956.mem", 138)
 Coverage.MallocInfo(16, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\initdt.jl.16956.mem", 156)
 Coverage.MallocInfo(32, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 463)
 Coverage.MallocInfo(48, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 214)
 Coverage.MallocInfo(48, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 215)
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 116)
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 195)
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 196)
 Coverage.MallocInfo(80, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 240)
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 287)
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 420)
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 438)
 Coverage.MallocInfo(96, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 457)
 Coverage.MallocInfo(128, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 224)
 Coverage.MallocInfo(144, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 295)
 Coverage.MallocInfo(144, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 390)
 Coverage.MallocInfo(208, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 213)
 Coverage.MallocInfo(224, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 267)
 Coverage.MallocInfo(448, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\perform_step\\low_order_rk_perform_step.jl.16956.mem", 564)
 Coverage.MallocInfo(496, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 262)
 Coverage.MallocInfo(832, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 66)
 Coverage.MallocInfo(1008, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\solve.jl.16956.mem", 339)
 Coverage.MallocInfo(3616, "C:\\Users\\sebastian\\.julia\\packages\\OrdinaryDiffEq\\UsWa2\\src\\perform_step\\low_order_rk_perform_step.jl.16956.mem", 599)
 ```

The largest allocation is in `perform_step` on line 599 for the error estimator
```
- @muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step=false)
0   @unpack t,dt,uprev,u,f,p = integrator
0   @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
0   k1 = integrator.fsalfirst
0   a = dt*a21
0   k2 = f(uprev+a*k1, p, t+c1*dt)
0   k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
0   k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
0   k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
0   g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
0   k6 = f(g6, p, t+dt)
0   u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
0   integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
0   integrator.destats.nf += 6
-   if typeof(integrator.alg) <: CompositeAlgorithm
-     g7 = u
-     # Hairer II, page 22
-     integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
-   end
0   if integrator.opts.adaptive
0     utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
0     atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
3616     integrator.EEst = integrator.opts.internalnorm(atmp,t)
-   end
0   integrator.k[1] = k1
0   integrator.k[2] = k2
0   integrator.k[3] = k3
0   integrator.k[4] = k4
0   integrator.k[5] = k5
0   integrator.k[6] = k6
0   integrator.k[7] = k7
0   integrator.u = u
- end
```
**Comparing with the previous case the allocation is about 2 times bigger.**

The second largest is in `__init` on line 339 for the integrator
```
- function DiffEqBase.__init(prob::Union{DiffEqBase.AbstractODEProblem,DiffEqBase.AbstractDAEProblem},
-                            alg::Union{OrdinaryDiffEqAlgorithm,DAEAlgorithm},
-                            timeseries_init = typeof(prob.u0)[],
-                            ts_init = eltype(prob.tspan)[],
-                            ks_init = [],
-                            recompile::Type{Val{recompile_flag}} = Val{true};
-                            saveat = eltype(prob.tspan)[],
-                            tstops = eltype(prob.tspan)[],
-                            d_discontinuities= eltype(prob.tspan)[],
-                            save_idxs = nothing,
-                            save_everystep = isempty(saveat),
-                            save_on = true,
-                            save_start = save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[1] in saveat,
-                            save_end = save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[2] in saveat,
-                            callback = nothing,
-                            dense = save_everystep && !(typeof(alg) <: Union{DAEAlgorithm,FunctionMap}) && isempty(saveat),
-                            calck = (callback !== nothing && callback != CallbackSet()) || # Empty callback
-                                    (!isempty(setdiff(saveat,tstops)) || dense), # and no dense output
-                            dt = typeof(alg) <: FunctionMap && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
-                            dtmin = typeof(one(eltype(prob.tspan))) <: AbstractFloat ? eps(eltype(prob.tspan)) :
-                                    typeof(one(eltype(prob.tspan))) <: Integer ? 0 : eltype(prob.tspan)(1//10^(10)),
-                            dtmax = eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
-                            force_dtmin = false,
-                            adaptive = isadaptive(alg),
-                            gamma = gamma_default(alg),
-                            abstol = nothing,
-                            reltol = nothing,
-                            qmin = qmin_default(alg),
-                            qmax = qmax_default(alg),
-                            qsteady_min = qsteady_min_default(alg),
-                            qsteady_max = qsteady_max_default(alg),
-                            qoldinit = 1//10^4,
-                            fullnormalize = true,
-                            failfactor = 2,
-                            beta1 = nothing,
-                            beta2 = nothing,
-                            maxiters = adaptive ? 1000000 : typemax(Int),
-                            internalnorm = ODE_DEFAULT_NORM,
-                            internalopnorm = LinearAlgebra.opnorm,
-                            isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
-                            unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
-                            verbose = true,
-                            timeseries_errors = true,
-                            dense_errors = false,
-                            advance_to_tstop = false,
-                            stop_at_next_tstop = false,
-                            initialize_save = true,
-                            progress = false,
-                            progress_steps = 1000,
-                            progress_name = "ODE",
-                            progress_message = ODE_DEFAULT_PROG_MESSAGE,
-                            userdata = nothing,
-                            allow_extrapolation = alg_extrapolates(alg),
-                            initialize_integrator = true,
-                            alias_u0 = false,
-                            kwargs...) where recompile_flag
-
832   if prob isa DiffEqBase.AbstractDAEProblem && alg isa OrdinaryDiffEqAlgorithm
-     error("You cannot use an ODE Algorithm with a DAEProblem")
-   end
-
-   if prob isa DiffEqBase.AbstractODEProblem && alg isa DAEAlgorithm
-     error("You cannot use an DAE Algorithm with a ODEProblem")
-   end
-
-   if typeof(prob.f)<:DynamicalODEFunction && typeof(prob.f.mass_matrix)<:Tuple
-     if any(mm != I for mm in prob.f.mass_matrix)
-       error("This solver is not able to use mass matrices.")
-     end
0   elseif !(typeof(prob)<:DiscreteProblem) &&
-          !(typeof(prob)<:DiffEqBase.AbstractDAEProblem) &&
-          !is_mass_matrix_alg(alg) &&
-          prob.f.mass_matrix != I
0     error("This solver is not able to use mass matrices.")
-   end
-
0   if !isempty(saveat) && dense
0     @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
-   end
-
0   progress && @logmsg(LogLevel(-1),progress_name,_id=_id = :OrdinaryDiffEq,progress=0)
-
0   tType = eltype(prob.tspan)
0   tspan = prob.tspan
0   tdir = sign(tspan[end]-tspan[1])
-
0   t = tspan[1]
-
0   if (((!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm) && !(typeof(alg) <: DAEAlgorithm)) || !adaptive) && dt == tType(0) && isempty(tstops)) && !(typeof(alg) <: Union{FunctionMap,LinearExponential})
0       error("Fixed timestep methods require a choice of dt or choosing the tstops")
-   end
-
0   f = prob.f
-   p = prob.p
-
-   # Get the control variables
-
0   if alias_u0
0     u = prob.u0
-   else
0     u = recursivecopy(prob.u0)
-   end
-
0   uType = typeof(u)
-   uBottomEltype = recursive_bottom_eltype(u)
-   uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)
-
80   ks = Vector{uType}(undef, 0)
-
-   uEltypeNoUnits = recursive_unitless_eltype(u)
-   tTypeNoUnits   = typeof(one(tType))
-
-   if typeof(alg) <: FunctionMap
-     abstol_internal = real.(zero.(u))
-   elseif abstol === nothing
-     if uBottomEltypeNoUnits == uBottomEltype
0       abstol_internal = real(convert(uBottomEltype,oneunit(uBottomEltype)*1//10^6))
-     else
-       abstol_internal = real.(oneunit.(u).*1//10^6)
-     end
-   else
-     abstol_internal = real.(abstol)
-   end
-
-   if typeof(alg) <: FunctionMap
-     reltol_internal = real.(zero(first(u)/t))
-   elseif reltol === nothing
-     if uBottomEltypeNoUnits == uBottomEltype
0       reltol_internal = real(convert(uBottomEltype,oneunit(uBottomEltype)*1//10^3))
-     else
-       reltol_internal = real.(oneunit.(u).*1//10^3)
-     end
-   else
-     reltol_internal = real.(reltol)
-   end
-
0   dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
-   # dtmin is all abs => does not care about sign already.
-
0   if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number && uBottomEltypeNoUnits == uBottomEltype # Could this be more efficient for other arrays?
-     if !(typeof(u) <: ArrayPartition)
-       rate_prototype = recursivecopy(u)
-     else
-       rate_prototype = similar(u, typeof.(oneunit.(recursive_bottom_eltype.(u.x))./oneunit(tType))...)
-     end
-   else
-     if uBottomEltypeNoUnits == uBottomEltype
-       rate_prototype = u
-     else # has units!
-       rate_prototype = u/oneunit(tType)
-     end
-   end
-   rateType = typeof(rate_prototype) ## Can be different if united
-
0   tstops_internal, saveat_internal, d_discontinuities_internal =
-     tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan)
-
-   callbacks_internal = CallbackSet(callback)
-
0   max_len_cb = DiffEqBase.max_vector_callback_length(callbacks_internal)
-   if max_len_cb isa VectorContinuousCallback
-     if isinplace(prob)
-       callback_cache = DiffEqBase.CallbackCache(u,max_len_cb.len,uBottomEltype,uBottomEltype)
-     else
-       callback_cache = DiffEqBase.CallbackCache(max_len_cb.len,uBottomEltype,uBottomEltype)
-     end
-   else
-     callback_cache = nothing
-   end
-
-   ### Algorithm-specific defaults ###
-   if save_idxs === nothing
-     ksEltype = Vector{rateType}
-   else
-     ks_prototype = rate_prototype[save_idxs]
-     ksEltype = Vector{typeof(ks_prototype)}
-   end
-
-   # Have to convert incase passed in wrong.
-   if save_idxs === nothing
-     timeseries = convert(Vector{uType},timeseries_init)
-   else
-     u_initial = u[save_idxs]
-     timeseries = convert(Vector{typeof(u_initial)},timeseries_init)
-   end
-   ts = convert(Vector{tType},ts_init)
80   ks = convert(Vector{ksEltype},ks_init)
80   alg_choice = Int[]
-
0   if !adaptive && save_everystep && tspan[2]-tspan[1] != Inf
0     dt == 0 ? steps = length(tstops) :
-               steps = ceil(Int,internalnorm((tspan[2]-tspan[1])/dt,tspan[1]))
0     sizehint!(timeseries,steps+1)
0     sizehint!(ts,steps+1)
0     sizehint!(ks,steps+1)
0   elseif save_everystep
0     sizehint!(timeseries,50)
0     sizehint!(ts,50)
0     sizehint!(ks,50)
0   elseif !isempty(saveat_internal)
0     sizehint!(timeseries,length(saveat_internal)+1)
0     sizehint!(ts,length(saveat_internal)+1)
0     sizehint!(ks,length(saveat_internal)+1)
-   else
208     sizehint!(timeseries,2)
48     sizehint!(ts,2)
48     sizehint!(ks,2)
-   end
-
0   if save_start
-     saveiter = 1 # Starts at 1 so first save is at 2
-     saveiter_dense = 1
0     copyat_or_push!(ts,1,t)
-     if save_idxs === nothing
0       copyat_or_push!(timeseries,1,u)
128       copyat_or_push!(ks,1,[rate_prototype])
-     else
-       copyat_or_push!(timeseries,1,u_initial,Val{false})
-       copyat_or_push!(ks,1,[ks_prototype])
-     end
-
-     if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
-       copyat_or_push!(alg_choice,1,1)
-     end
-   else
-     saveiter = 0 # Starts at 0 so first save is at 1
-     saveiter_dense = 0
-   end
-
0   QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits
-
80   k = rateType[]
-
-   if uses_uprev(alg, adaptive) || calck
0     uprev = recursivecopy(u)
-   else
-     # Some algorithms do not use `uprev` explicitly. In that case, we can save
-     # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
-     uprev = u
-   end
0   if allow_extrapolation
0     uprev2 = recursivecopy(u)
-   else
-     uprev2 = uprev
-   end
-
0   cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol_internal,p,calck,Val(isinplace(prob)))
-
-   if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
-     id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,dense,cache)
-     beta2 === nothing && ( beta2=beta2_default(alg.algs[cache.current]) )
-     beta1 === nothing && ( beta1=beta1_default(alg.algs[cache.current],beta2) )
-   else
496     id = InterpolationData(f,timeseries,ts,ks,dense,cache)
0     beta2 === nothing && ( beta2=beta2_default(alg) )
0     beta1 === nothing && ( beta1=beta1_default(alg,beta2) )
-   end
-
224   opts = DEOptions{typeof(abstol_internal),typeof(reltol_internal),QT,tType,
-                    typeof(internalnorm),typeof(internalopnorm),typeof(callbacks_internal),typeof(isoutofdomain),
-                    typeof(progress_message),typeof(unstable_check),typeof(tstops_internal),
-                    typeof(d_discontinuities_internal),typeof(userdata),typeof(save_idxs),
-                    typeof(maxiters),typeof(tstops),typeof(saveat),
-                    typeof(d_discontinuities)}(
-                        maxiters,save_everystep,adaptive,abstol_internal,
-                        reltol_internal,QT(gamma),QT(qmax),
-                        QT(qmin),QT(qsteady_max),
-                        QT(qsteady_min),QT(failfactor),tType(dtmax),
-                        tType(dtmin),internalnorm,internalopnorm,save_idxs,tstops_internal,saveat_internal,
-                        d_discontinuities_internal,
-                        tstops,saveat,d_discontinuities,
-                        userdata,progress,progress_steps,
-                        progress_name,progress_message,timeseries_errors,dense_errors,
-                        QT(beta1),QT(beta2),QT(qoldinit),dense,
-                        save_on,save_start,save_end,callbacks_internal,isoutofdomain,
-                        unstable_check,verbose,
-                        calck,force_dtmin,advance_to_tstop,stop_at_next_tstop)
-
96   destats = DiffEqBase.DEStats(0)
-
-   if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
-     sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,
-                       dense=dense,k=ks,interp=id,
-                       alg_choice=alg_choice,
-                       calculate_error = false, destats=destats)
-   else
144     sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,
-                       dense=dense,k=ks,interp=id,
-                       calculate_error = false, destats=destats)
-   end
-
-   if recompile_flag == true
-     FType = typeof(f)
-     SolType = typeof(sol)
-     cacheType = typeof(cache)
-   else
-     FType = Function
-     if alg isa OrdinaryDiffEqAlgorithm
-       SolType = DiffEqBase.AbstractODESolution
-       cacheType =  OrdinaryDiffEqCache
-     else
-       SolType = DiffEqBase.AbstractDAESolution
-       cacheType =  DAECache
-     end
-   end
-
-   # rate/state = (state/time)/state = 1/t units, internalnorm drops units
-   eigen_est = one(uBottomEltypeNoUnits)/one(tType)
-   tprev = t
-   dtcache = tType(dt)
-   dtpropose = tType(dt)
-   iter = 0
-   kshortsize = 0
-   reeval_fsal = false
-   u_modified = false
-   EEst = tTypeNoUnits(1)
-   just_hit_tstop = false
-   isout = false
-   accept_step = false
-   force_stepfail = false
-   last_stepfail = false
-   event_last_time = 0
-   vector_event_last_time = 1
-   last_event_error = zero(uBottomEltypeNoUnits)
-   dtchangeable = isdtchangeable(alg)
-   q11 = tTypeNoUnits(1)
-   success_iter = 0
-   erracc = tTypeNoUnits(1)
-   dtacc = tType(1)
-
1008   integrator = ODEIntegrator{typeof(alg),isinplace(prob),uType,tType,typeof(p),typeof(eigen_est),
-                              QT,typeof(tdir),typeof(k),SolType,
-                              FType,cacheType,
-                              typeof(opts),fsal_typeof(alg,rate_prototype),
-                              typeof(last_event_error),typeof(callback_cache)}(
-                              sol,u,k,t,tType(dt),f,p,uprev,uprev2,tprev,
-                              alg,dtcache,dtchangeable,
-                              dtpropose,tdir,eigen_est,EEst,QT(qoldinit),q11,
-                              erracc,dtacc,success_iter,
-                              iter,saveiter,saveiter_dense,cache,callback_cache,
-                              kshortsize,force_stepfail,last_stepfail,
-                              just_hit_tstop,event_last_time,vector_event_last_time,last_event_error,
-                              accept_step,
-                              isout,reeval_fsal,
-                              u_modified,opts,destats)
0   if initialize_integrator
0     initialize_callbacks!(integrator, initialize_save)
0     initialize!(integrator,integrator.cache)
0     save_start && typeof(alg) <: CompositeAlgorithm && copyat_or_push!(alg_choice,1,integrator.cache.current)
-   end
-
0   handle_dt!(integrator)
-
0   integrator
- end
```

Next follows line 66 in the same function (see above), and I think it is associated
with the creation of `prob` and next is 262 with interpolation data(`id`).

Next follows on line 564 in the `initialize!` function
```
- function initialize!(integrator, cache::Tsit5ConstantCache)
0   integrator.kshortsize = 7
448   integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
0   integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
0   integrator.destats.nf += 1
-
-   # Avoid undefined entries if k is an array of arrays
0   integrator.fsallast = zero(integrator.fsalfirst)
0   integrator.k[1] = integrator.fsalfirst
0   @inbounds for i in 2:integrator.kshortsize-1
0     integrator.k[i] = zero(integrator.fsalfirst)
-   end
0   integrator.k[integrator.kshortsize] = integrator.fsallast
- end
```
As we can observe, the allocations remained the same for the rest.

Running the last case and disabling inlining gives
```
& "C:\Program Files\Julia-1.3.0-rc3\bin\julia.exe" --track-allocation=user --inline=no .\test.jl
  9.043187 seconds (7.23 M allocations: 390.388 MiB, 7.61% gc time)
  1.343243 seconds (2.04 M allocations: 99.522 MiB, 0.54% gc time)
  1.406408 seconds (2.04 M allocations: 99.522 MiB, 0.73% gc time)
elapsed time (ns): 1406407900
gc time (ns):      10232400
bytes allocated:   104356704
pool allocs:       2044948
GC pauses:         3
```
Comparing with the very first case, we can se a roughly 2x increase in the allocations.
In the user function we have
```
        - function lorenz2(u,p,t)
    86976     du1 = 10.0*(u[2]-u[1])
   130464     du2 = u[1]*(28.0-u[3]) - u[2]
   152208     du3 = u[1]*u[2] - (8/3)*u[3]
        -
    86976     du4 = 10.0*(u[2+3]-u[1+3])
   130464     du5 = u[1+3]*(28.0-u[3+3]) - u[2+3]
   152208     du6 = u[1+3]*u[2+3] - (8/3)*u[3+3]
        -
        0     return SVector{6}(du1,du2,du3,du4,du5,du6)
        - end
```
in solve
```
        - function DiffEqBase.solve!(integrator::ODEIntegrator)
        0   @inbounds while !isempty(integrator.opts.tstops)
    14528     while integrator.tdir * integrator.t < top(integrator.opts.tstops)
        0       loopheader!(integrator)
        0       if check_error!(integrator) != :Success
        0         return integrator.sol
        -       end
   126560       perform_step!(integrator,integrator.cache)
        0       loopfooter!(integrator)
        0       if isempty(integrator.opts.tstops)
        -         break
        -       end
        -     end
        0     handle_tstop!(integrator)
        -   end
        0   postamble!(integrator)
        -
       96   f = integrator.sol.prob.f
        -
        0   if DiffEqBase.has_analytic(f)
        -     DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
        -   end
        0   if integrator.sol.retcode != :Default
        0     return integrator.sol
        -   end
        0   integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
        - end
```
and perform_step
```
        - @muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step=false)
   166336   @unpack t,dt,uprev,u,f,p = integrator
   122944   @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
    14464   k1 = integrator.fsalfirst
     3616   a = dt*a21
    32544   k2 = f(uprev+a*k1, p, t+c1*dt)
    75936   k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
   104864   k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
   133792   k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
   159104   g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
     3616   k6 = f(g6, p, t+dt)
   188032   u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    32544   integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
     9072   integrator.destats.nf += 6
        0   if typeof(integrator.alg) <: CompositeAlgorithm
        -     g7 = u
        -     # Hairer II, page 22
        -     integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
        -   end
        0   if integrator.opts.adaptive
   202496     utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
     7232     atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    10848     integrator.EEst = integrator.opts.internalnorm(atmp,t)
        -   end
    14464   integrator.k[1] = k1
    14464   integrator.k[2] = k2
    14464   integrator.k[3] = k3
    14464   integrator.k[4] = k4
    14464   integrator.k[5] = k5
    14464   integrator.k[6] = k6
    14464   integrator.k[7] = k7
    14464   integrator.u = u
        - end
```

Again comparing, there are about 2x more allocations.
