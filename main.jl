using Plots, OrdinaryDiffEq, StaticArrays, Statistics, BlackBoxOptim

# Get data for trial
include("getdata.jl")
fname = "MPIData/MPI-Advanced-Ethan0030.mat"
splx, sply, splz, base, mid, tip, initconds, la, lb, tspan = getdata(fname)

# Create model
include("parameters.jl")
include("initconds.jl")
include("model/eom.jl")

pfixed = (splx, sply, splz, la, lb)
pvariable = (ka = 10.0, kb = 10.0, ba = 0.1, bb = 0.1, ma = 0.5, mb = 0.5)
p = getparameters((pfixed..., pvariable...))
u₀ = getinitcond(initconds, p)

prob = ODEProblem(eom!, u₀, tspan, p)
# sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-7)

# Find optimal parameters
include("cost.jl")
include("model/functions.jl")
times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
bounds = [(0, 10), (0, 10), (0, 10), (0, 10), (0, 5), (0, 5)]

res = bboptimize(cost, SearchRange=bounds, NumDimensions=6) # , NThreads=Threads.nthreads()-1)
opt = best_candidate(res)
score = best_fitness(res)

# Remake problem with optimal parameters
pvariable = (ka = opt[1], kb = opt[2], ba = opt[3], bb = opt[4], ma = opt[5], mb = opt[6])
pnew = merge(p, pvariable)
# newprob = remake(prob, p=pnew)
prob = ODEProblem(eom!, u₀, tspan, pnew)


# Solve
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-7)

# Get data for comparison
base_sim = [splx(times) sply(times) splz(times)]
mid_sim = [tup[j] for tup ∈ [p2(sol, t) for t ∈ times], j ∈ 1:3]
tip_sim = [tup[j] for tup ∈ [p3(sol, t) for t ∈ times], j ∈ 1:3]

# Save to file
name = split(fname, "-") |> x -> split(x[3], ".mat")[1]
open("result.csv", "a") do io
    write(io, name, opt, score)
end


# Plot
include("plotting.jl")
plot_comparison(mid,tip,mid_sim,tip_sim)
animate_comparison(base,mid,tip,base_sim,mid_sim,tip_sim)
