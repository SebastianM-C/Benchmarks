using StaticArrays
using Profile
using BenchmarkTools
using OrdinaryDiffEq

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

@inbounds @inline function ṗ(p, q, params, t)
    A, B, D = params
    dp1 = -A * q[1] - 3 * B / √2 * (q[2]^2 - q[1]^2) - D * q[1] * (q[1]^2 + q[2]^2)
    dp2 = -q[2] * (A + 3 * √2 * B * q[1] + D * (q[1]^2 + q[2]^2))
    return SVector{2}(dp1, dp2)
end

@inline function q̇(p, q, params, t)
    params.A * p
end

const q0 = SVector{2}([0.0, -4.363920590485035])
const p0 = SVector{2}([10.923918825236079, -5.393598858645495])
const z0 = vcat(p0, q0)
const p = (A=1,B=0.55,D=0.4)

const tspan = (0., 10.)

prob1 = ODEProblem(ż, z0, tspan, p)

solve(prob1, Vern9(), abstol=1e-14, reltol=1e-14);
@timev solve(prob1, Vern9(), abstol=1e-14, reltol=1e-14);
Profile.clear_malloc_data()
solve(prob1, Vern9(), abstol=1e-14, reltol=1e-14);
exit()

# ṗ(p0, q0, p, 1.);
# Profile.clear_malloc_data()
# @timev ṗ(p0, q0, p, 1.);
#
# # @benchmark q̇(p0, q0, p, 1.)
#
# @timev ż(z0, p, 1.);
# Profile.clear_malloc_data()
# @timev ż(z0, p, 1.);
