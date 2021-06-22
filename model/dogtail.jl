using OrdinaryDiffEq, Plots

include("eom.jl")
include("functions.jl")

# Parameters
ma = 1; mb = 1;
la = lb = 0.5; lao = lbo = 0.25
ixa = ixb = 0.01;
iya = iza = 1 / 12 * ma * lao; iyb = izb = 1 / 12 * mb * lbo
g = -9.81
ka = 1; kb = 10; ba = bb = 0.1

# Specified base function 
pf(t)   = 0, 0, 0 # 0.1 * sin(10 * t)
pfp(t)  = 0, 0, 0 # cos(10 * t)
pfpp(t) = 0, 0, 0 # -10 * sin(10 * t)

p = (ma, ixa, iya, iza, lao, la, mb, ixb, iyb, izb, lbo, lb, g, ka, kb, ba, bb, pf = pf, pfp = pfp, pfpp = pfpp)


# Initial conditions
q1₀ = 0; q2₀ = 0; q3₀ = 0; q4₀ = 0; q5₀ = 0; q6₀ = 0
u1₀ = 0; u2₀ = 100; u3₀ = 100; u4₀ = 0; u5₀ = 100; u6₀ = 0

u₀ = deg2rad.([q1₀, q2₀, q3₀, q4₀, q5₀, q6₀, u1₀, u2₀, u3₀, u4₀, u5₀, u6₀])

# Time span
t₀ = 0.0
t₁ = 5.0
tspan = (t₀, t₁)

# Set up problem
prob = ODEProblem(eom!, u₀, tspan, p)

# Solve
sol = solve(prob, Tsit5(), saveat=0.01, abstol=1e-8, reltol=1e-8)

# Plot and animate
plot_pendulum(sol, 0, tip=true)
plot_energy(sol)
animate_pendulum(sol)