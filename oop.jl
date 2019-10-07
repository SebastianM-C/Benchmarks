using OrdinaryDiffEq
using BenchmarkTools
using StaticArrays

function lorenz(u,p,t)
  du1 = 10.0*(u[2]-u[1])
  du2 = u[1]*(28.0-u[3]) - u[2]
  du3 = u[1]*u[2] - (8/3)*u[3]

  return SVector{3}(du1,du2,du3)
end

@inbounds @inline function ż(z, p, t)
    A, B, D = p
    p₀, p₂ = z[SVector{2}(1:2)]
    q₀, q₂ = z[SVector{2}(3:4)]

    return SVector{4}(
        -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2),
        -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2)),
        A * p₀,
        A * p₂
    )
end

henon(z, p, t) = SVector(
        -z[3] * (1 + 2z[4]),
        -z[4] - (z[3]^2 - z[4]^2),
        z[1],
        z[2]
    )

u0 = @SVector [1.0,0.0,0.0]
u = vcat(u0,u0)
p = (A=1,B=0.55,D=0.4)

q0 = SVector{2}([0.0, -4])
p0 = SVector{2}([10, -5])
z0 = vcat(p0, q0)
z1 = SVector{4}([0.1,0,0,0.5])

tspan1 = (0.0,10.0)
prob1_ok = ODEProblem(lorenz,u0,tspan1)
prob1_notok = ODEProblem(ż,z0,tspan1,p)
prob1 = ODEProblem(henon,z1,tspan1)
@btime solve($prob1_ok, Vern9(), save_everystep=false)
@btime solve($prob1_notok, Vern9(), save_everystep=false)
@btime solve($prob1, Vern9(), save_everystep=false)

tspan2 = (0.0,100.0)
prob2_ok = ODEProblem(lorenz,u0,tspan2)
prob2_notok = ODEProblem(ż,z0,tspan2,p)
prob2 = ODEProblem(henon,z1,tspan2)

@btime solve($prob2_ok, Vern9(), save_everystep=false)
@btime solve($prob2_notok, Vern9(), save_everystep=false)
@btime solve($prob2, Vern9(), save_everystep=false)

@btime ż($z0,$p,0)
@code_warntype ż(z0,p,0)
