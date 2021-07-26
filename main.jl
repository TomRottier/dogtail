using Plots, OrdinaryDiffEq, StaticArrays, BlackBoxOptim, Parameters, Evolutionary, SimulatedAnnealing
using DelimitedFiles, Statistics, LinearAlgebra, Distributed

# Parallel version
# using Plots
# using Distributed
# addprocs(4)
# @everywhere begin
#     # using Pkg; Pkg.activate(".")
#     using OrdinaryDiffEq, StaticArrays, Parameters, SimulatedAnnealing
#     using DelimitedFiles, StaticArrays, Statistics, LinearAlgebra

#     include("getdata.jl")
#     include("parameters.jl")
#     include("orientations.jl")
#     include("initconds.jl")
#     include("model/eom.jl")
#     include("cost.jl")
#     include("model/functions.jl")
# end

include("getdata.jl")
include("parameters.jl")
include("orientations.jl")
include("initconds.jl")
include("model/eom.jl")
include("cost.jl")
include("model/functions.jl")
include("plotting.jl")

fnames = readdir("Data/")
mids_sim = Vector{Matrix{Float64}}(undef, length(fnames))
tips_sim = similar(mids_sim)
mids = similar(mids_sim)
tips = similar(mids)
dognames = Vector{String}(undef, length(fnames))

#= Threads.@threads =# for (i, fname) in collect(enumerate(fnames))
    # Get data for trial
    println(i)
    @everywhere splx, sply, splz, time, base, mid, tip, orientations, la, lb, tspan = getdata("Data/" * $fname)
    mids[i] = mid; tips[i] = tip
    dognames[i] = split(fname, ".")[1]
    
    # Create model
    @everywhere begin
        p = initialise_parameters(splx, sply, splz, la, lb)
        times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
        # u₀ = [orientations[i][1] for i ∈ eachindex(orientations)]
        # prob = ODEProblem(eom!, u₀, tspan, p)
        orientations[10] .= 0.0
        u₀_SA = SVector{12,Float64}([orientations[j][1] for j ∈ eachindex(orientations)])
        prob = ODEProblem(eom_SA, u₀_SA, tspan, p)        
    end

    # Find optimal parameters

    # Using BlackBoxOptim
    bounds = [(0.05, 0.6), (0.05, 0.6), (0.0001, 0.01), (0.0002, 0.01), (1e-6, 0.01), (1e-6, 0.01), (-π, π), (-π, π), (-π, π)]
    res = bboptimize(x -> cost(x, p, prob, times, mid, tip), SearchRange=bounds, NumDimensions=length(bounds), MaxFuncEvals=100, ftol=1e-4)
    opt = best_candidate(res)
    score = best_fitness(res)

    # Using Evolutionary.jl
    # lb = [0.05, 0.05, 0.0001, 0.0001, 1e-6, 1e-6, -π, -π, -π]
    # ub = [0.5, 0.5, 0.005, 0.005, 0.0002, 0.0002, π, π, π]
    # x₀ = @. lb + (ub - lb) / 2
    # ga = GA(populationSize=100)
    # opts = Evolutionary.Options(show_trace=true)
    # res = Evolutionary.optimize(x -> cost(x, p, prob, time, mid, tip), lb, ub, x₀, ga, opts)

    # Using my SimulatedAnnealing
    # lb = [0.05, 0.05, 0.0001, 0.0001, 1e-6, 1e-6, -π, -π, -π]
    # ub = [0.5, 0.5, 0.005, 0.005, 0.0002, 0.0002, π, π, π]
    # x₀ = @. lb + (ub - lb) / 2
    # options = Options(func=x -> cost(x, p, prob, times, mid, tip), N=9, lb=lb, ub=ub, tol=1e-3, print_status=true)
    # current = State(f=options.f(x₀), x=x₀, v=ub .- lb, T=0.001)
    # result = Result(options.f(x₀), x₀)
    # sa!(current, result, options)


    # Call cost with optimal paramters
    mid_sim, tip_sim = cost(opt, p, prob, times)
    mids_sim[i] = mid_sim; tips_sim[i] = tip_sim

    # Save to file
    # name = split(fname, ".")[1]
    # open("results.csv", "a") do io
    #     writedlm(io, [name  opt... score], ',')
    # end
    # open("SimData/" * name * ".csv", "w") do io
    #     header = ["time" "midX" "midY" "midZ" "tipX" "tipY" "tipZ"]
    #     writedlm(io, [ header; times mid_sim tip_sim], ',')
    # end

end

# Plot
for i ∈ eachindex(mids)
    mid = mids[i]
    tip = tips[i]
    mid_sim = mids_sim[i]
    tip_sim = tips_sim[i]

    title = plot(title=dognames[i], grid=false, showaxis=false, ticks=false, bottom_margin=-50Plots.px)
    plt = plot_comparison(mid, tip, mid_sim, tip_sim)
    bigplt = plot(title, plt, layout=@layout([A{0.15h}; B]))
    display(bigplt)
end


results = readdlm("results.csv", ',', skipstart=13)

for result ∈ eachrow(results)
    name = result[1] * ".csv"
    opt = result[2:end]
    splx, sply, splz, time, base, mid, tip, orientations, la, lb, tspan = getdata("Data/" * name)

    p = initialise_parameters(splx, sply, splz, la, lb)
    times = collect(range(tspan[1], tspan[2], length=size(base, 1)))
    u₀_SA = SVector{12,Float64}([orientations[j][1] for j ∈ eachindex(orientations)])
    prob = ODEProblem(eom_SA, u₀_SA, tspan, p)
    mid_sim, tip_sim = cost(opt, p, prob, times)

    open("SimData/" * name, "w") do io
        header = ["time" "midX" "midY" "midZ" "tipX" "tipY" "tipZ"]
        writedlm(io, [ header; times mid_sim tip_sim], ',')
    end
end