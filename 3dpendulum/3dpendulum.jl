using Plots, OrdinaryDiffEq

include("eom.jl")
include("functions.jl")

# Parameters
m = 10
ix = 1; iy = 1; iz = 1;
lao = 0.5; la = 1.0
g = -9.81

p = (m, ix, iy, iz, lao, la, g)

# Initial conditions
q4₀ = 0; q5₀ = 0; q6₀ = 90;
u4₀ = 0; u5₀ = 0; u6₀ = 0;
u₀ = deg2rad.([q4₀, q5₀, q6₀, u4₀, u5₀, u6₀])

# Time span
tinitial = 0.0;
tfinal = 5.0;
tspan = (tinitial, tfinal)

# Set up problem
prob = ODEProblem(eom, u₀, tspan, p)

# Solve
sol = solve(prob, Tsit5(), saveat=0.001, abstol=1e-8, reltol=1e-8)


# plot_pendulum(sol, 0)
animate_pendulum(sol)
plot(totalenergy(sol))
plot!(kineticenergy(sol))
plot!(potentialenergy(sol))
