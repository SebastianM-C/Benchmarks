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
```
const tspan = (0.0,10.0)
const prob = ODEProblem(lorenz,u0,tspan)
# const prob = ODEProblem(lorenz2,u,tspan)
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
Strangely no allocations appear in `perform_step!`
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
let us consider
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
which has roughly double the allocations, but again the allocations are not shown in `perform_step`.
This may indicate that the problem still exists, but it is more difficult to debug due to inlining.

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
