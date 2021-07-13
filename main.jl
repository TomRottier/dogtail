using Plots, OrdinaryDiffEq, StaticArrays, BlackBoxOptim, Parameters, Evolutionary
using DelimitedFiles, Statistics, LinearAlgebra

include("getdata.jl")
include("parameters.jl")
include("orientations.jl")
include("initconds.jl")
include("model/eom.jl")
include("cost.jl")
include("model/functions.jl")
include("plotting.jl")

fnames = readdir("MPIData/")
sols = Vector(undef, length(fnames))
data = similar(sols)
plts = similar(sols)
names = similar(sols)

#= Threads.@threads =# for (i, fname) in collect(enumerate(fnames))
    # Get data for trial
    # fname = "MPIData/MPI-Advanced-Ethan0030.mat"
    println(i)
    splx, sply, splz, base, mid, tip, orientations, la, lb, tspan = getdata("MPIData/" * fname)
    # data[i] = cat(base, mid, tip, dims=3)
    # times = collect(range(tspan[1], tspan[2], length=size(base, 1)))

    # plt = plot(derivative(splx, times, 2))
    # plot!(derivative(sply, times, 2))
    # plot!(derivative(splz, times, 2), label=["x" "y" "z"], title="$i")
    # display(plt)

    # plt2 = plot(evaluate(splx, times), evaluate(splx, times), evaluate(splx, times), label="", title="$i")
    # display(plt2)

    # Create model
    p = initialise_parameters(splx, sply, splz, la, lb)
    u₀ = [orientations[i][1] for i ∈ eachindex(orientations)]
    times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
    prob = ODEProblem(eom!, u₀, tspan, p)

    # Find optimal parameters

    # Using BlackBoxOptim
    # bounds = [(0.05, 0.5), (0.05, 0.5), (0.0001, 0.005), (0.0001, 0.005), (1e-6, 0.0002), (1e-6, 0.0002), (-π, π), (-π, π), (-π, π)]
    # res = bboptimize(x -> cost(x, p, prob, times, mid, tip), SearchRange=bounds, NumDimensions=length(bounds), MaxFuncEvals=10)
    # opt = best_candidate(res)
    # score = best_fitness(res)

    # Using Evolutionary.jl
    lb = [0.05, 0.05, 0.0001, 0.0001, 1e-6, 1e-6, -π, -π, -π]
    ub = [0.5, 0.5, 0.005, 0.005, 0.0002, 0.0002, π, π, π]
    x₀ = @. lb + (ub - lb) / 2
    ga = GA(populationSize=100, selection=uniformranking(3), mutation=gaussian(), crossover=uniformbin())
    opts = Evolutionary.Options(show_trace=true)
    res = Evolutionary.optimize(x -> cost(x, p, prob, time, mid, tip), lb, ub, x₀, ga, opts)

    # Call cost with optimal paramters
    mid_sim, tip_sim = cost(opt, p, u₀, tspan, times)

    # # Save to file
    # name = split(fname, "-") |> x -> split(x[3], ".mat")[1]
    # names[i] = name
    # open("results.csv", "a") do io
    #     writedlm(io, [name  opt... score], ',')
    # end

end

# # Plot
for i ∈ eachindex(sols)
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


    plt = plot_comparison(mid, tip, mid_sim, tip_sim)
    plot!(title=names[i])
    display(plt)
end
