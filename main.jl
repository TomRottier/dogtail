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
coef = rand(MMatrix{6,6,Float64})
rhs = rand(MVector{6,Float64})

# prob = ODEProblem(eom!, u₀, tspan, p)

# Find optimal parameters
include("cost.jl")
include("model/functions.jl")
times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
bounds = [(0, 10), (0, 10), (0, 10), (0, 10), (0, 5), (0, 5)]
res = bboptimize(cost, SearchRange=bounds, NumDimensions=6)


# x₀ = [i for i in pvariable]
# options = Options(func=cost, N=6, Ns=20, Nt=5, lb=repeat([0.], 6), ub=repeat([10.], 6), print_status=true, c=repeat([2.0], 6))
# result = Result(fopt=options.f(x₀), xopt=x₀)
# current = State(f=options.f(x₀), x=x₀, v=options.ub .- options.lb, T=100.0)

# sa!(current, result, options)

opt = result.xopt

# Remake problem with optimal parameters
pvariable = (ka = pin[1], kb = pin[2], ba = pin[3], bb = pin[4], ma = pin[5], mb = pin[6])
pnew = merge(p, pvariable)
# newprob = remake(prob, p=pnew)
prob = ODEProblem(eom!, u₀, tspan, pnew)


# Solve
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-7)

# Get data for comparison
base_sim = [splx(times) sply(times) splz(times)]
mid_sim = [tup[j] for tup ∈ [p2(sol, t) for t ∈ times], j ∈ 1:3]
tip_sim = [tup[j] for tup ∈ [p3(sol, t) for t ∈ times], j ∈ 1:3]


# Plot
include("plotting.jl")
plot_comparison(mid,tip,mid_sim,tip_sim)
animate_comparison(base,mid,tip,base_sim,mid_sim,tip_sim)
