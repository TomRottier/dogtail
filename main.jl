using Plots, OrdinaryDiffEq, StaticArrays, Statistics, BlackBoxOptim, DelimitedFiles, Dierckx, MAT

include("getdata.jl")
include("parameters.jl")
include("initconds.jl")
include("model/eom.jl")
include("cost.jl")
include("model/functions.jl")
include("plotting.jl")

fnames = readdir("MPIData/")
sols = Vector(undef, length(fnames))
data = similar(sols)
plts = similar(sols)
names= similar(sols)

Threads.@threads for (i,fname) in collect(enumerate(fnames))
    # Get data for trial
    # fname = "MPIData/MPI-Advanced-Ethan0030.mat"
    splx, sply, splz, base, mid, tip, initconds, la, lb, tspan = getdata("MPIData/"*fname)
    data[i] = cat(base,mid,tip, dims=3)

    # Create model
    pfixed = (splx, sply, splz, la, lb)
    pvariable = (10,10,0.1,0.1,5,5)
    p = getparameters((pfixed..., pvariable...))
    u₀ = try 
        getinitcond(initconds, p)
    catch e
        println(e)
        continue
    end
    prob = ODEProblem(eom!, u₀, tspan, p)

    # Find optimal parameters
    times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
    bounds = [(0, 50), (0, 50), (0, 50), (0, 50), (0, 2), (0, 2)]
    res = bboptimize(x->cost(x, p, u₀, tspan, times, mid, tip), SearchRange=bounds, NumDimensions=6, MaxFuncEvals=1000)
    opt = best_candidate(res)
    score = best_fitness(res)

    # Remake problem with optimal parameters
    pnew = (p[1:end-5]..., opt...)
    prob = ODEProblem(eom!, u₀, tspan, pnew)

    # Solve
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-7)
    sols[i] = sol

    # Get data for comparison
    base_sim = [splx(times) sply(times) splz(times)]
    mid_sim = [tup[j] for tup ∈ [p2(sol, t) for t ∈ times], j ∈ 1:3]
    tip_sim = [tup[j] for tup ∈ [p3(sol, t) for t ∈ times], j ∈ 1:3]

    # Save to file
    name = split(fname, "-") |> x -> split(x[3], ".mat")[1]
    names[i] = name
    open("results.csv", "a") do io
        writedlm(io, [name  opt... score], ',')
    end

end

# Plot
for  i in 1:length(sols)
    sol = try
        sols[i]
    catch
        continue
    end

    @inbounds la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb = sol.prob.p
    tspan = sol.prob.tspan

    base = data[i][:,:,1]
    mid = data[i][:,:,2]
    tip = data[i][:,:,3]

    times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
    base_sim = [fx(times) fy(times) fz(times)]
    mid_sim = [tup[j] for tup ∈ [p2(sol, t) for t ∈ times], j ∈ 1:3]
    tip_sim = [tup[j] for tup ∈ [p3(sol, t) for t ∈ times], j ∈ 1:3]


    plt = plot_comparison(mid,tip,mid_sim,tip_sim)
    plot!(title=names[i])
    display(plt)
end
