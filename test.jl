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
# const prob = ODEProblem(lorenz,u0,tspan)
const prob = ODEProblem(lorenz2,u,tspan)

@time solve(prob, Tsit5(), save_everystep=false);

@time solve(prob, Tsit5(), save_everystep=false);
Profile.clear_malloc_data()
@timev solve(prob, Tsit5(), save_everystep=false);
exit()
