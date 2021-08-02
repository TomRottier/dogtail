using Plots, OrdinaryDiffEq, StaticArrays, BlackBoxOptim, Parameters, SimulatedAnnealing, Dierckx
using DelimitedFiles, Statistics, LinearAlgebra, Distributed

include("model/createModel.jl")
include("orientations.jl")
include("getdata.jl")
include("parameters.jl")
include("cost.jl")
include("model/functions.jl")
include("plotting.jl")

fnames = readdir("Data/")
mids_sim = Vector{Matrix{Float64}}(undef, length(fnames))
tips_sim = similar(mids_sim)
mids = similar(mids_sim)
tips = similar(mids)
dognames = Vector{String}(undef, length(fnames))

Threads.@threads for (i, fname) in collect(enumerate(fnames))
    # Get data for trial and create model
    println(i)
    prob, p, time, base, mid, tip = createModel("Data//" * fname)
    mids[i] = mid; tips[i] = tip
    dognames[i] = split(fname, ".")[1]
    
    # Find optimal parameters

    # BlackBoxOptim
    # bounds = [(0.01, 0.6), (0.01, 0.6), (0.0001, 0.1), (0.0002, 0.1), (1e-6, 0.1), (1e-6, 0.1), (-π, π), (-π, π), (-π, π)]
    # res = bboptimize(x -> cost(x, p, prob, time, mid, tip), SearchRange=bounds, NumDimensions=length(bounds), MaxFuncEvals=1000, ftol=1e-4)
    # opt = best_candidate(res)
    # score = best_fitness(res)

    # SimulatedAnnealing
    lb = [0.01, 0.01, 0.0001, 0.0001, 1e-6, 1e-6, -π, -π, -π]
    ub = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, π, π, π]
    x₀ = @. lb + (ub - lb) / 2
    options = Options(func=x -> cost(x, p, prob, time, mid, tip), N=9, Nt=5, lb=lb, ub=ub, tol=1e-3, print_status=true)
    current = State(f=options.f(x₀), x=x₀, v=ub .- lb, T=0.0005)
    result = Result(options.f(x₀), x₀)
    sa!(current, result, options)
    opt = result.xopt
    score = result.fopt


    # Call cost with optimal paramters
    mid_sim, tip_sim = cost(opt, p, prob, time)
    mids_sim[i] = mid_sim; tips_sim[i] = tip_sim

    # Save to file
    name = split(fname, ".")[1]
    open("results.csv", "a") do io
        writedlm(io, [name  opt... score], ',')
    end
    open("SimData/" * name * ".csv", "w") do io
        header = ["time" "midX" "midY" "midZ" "tipX" "tipY" "tipZ"]
        writedlm(io, [ header; time mid_sim tip_sim], ',')
    end

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
