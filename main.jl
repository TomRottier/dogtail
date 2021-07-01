using Plots, OrdinaryDiffEq, StaticArrays, Statistics

# Get data for trial
include("getdata.jl")
fname = "MPIData/MPI-Advanced-Ethan0030.mat"

splx, sply, splz, base, mid, tip, initconds, la, lb, tspan = getdata(fname)

# Create simulation and run
include("parameters.jl")
include("initconds.jl")
include("model/eom.jl")

p = parameters(;fx=splx, fy=sply, fz=splz, ka=0.1, kb=0.01, ba=0.1, bb=0.1, la=la, lb=lb)
u₀ = getinitcond(initconds, p)
coef = rand(MMatrix{6,6,Float64})
rhs = rand(MVector{6,Float64})

prob = ODEProblem(eom!, u₀, tspan, p)
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-7)

# Compare simulation output to measured data
include("model/functions.jl")
time = collect(range(tspan[1], tspan[2], length=size(base, 1)))
base_sim = [splx(time) sply(time) splz(time)]
mid_sim = [tup[j] for tup ∈ [p2(sol, t) for t ∈ time], j ∈ 1:3]
tip_sim = [tup[j] for tup ∈ [p3(sol, t) for t ∈ time], j ∈ 1:3]

# Average RMSE from each direction
mid_rmse = sqrt.((mean((mid - mid_sim).^2, dims=1))) |> mean
tip_rmse = sqrt.((mean((tip - tip_sim).^2, dims=1))) |> mean


# Plot
include("plotting.jl")
plot_comparison(mid,tip,mid_sim,tip_sim)
animate_comparison(base,mid,tip,base_sim,mid_sim,tip_sim)
