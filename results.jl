using Plots
using DelimitedFiles, Statistics


include("plotting.jl")

# Load data and plot
fnames = readdir("SimData/")
plts = Vector(undef, length(fnames))

for (idx, fname) ∈ enumerate(fnames)
    fname_sim = "SimData/" * fname
    fname_exp = "Data/" * fname
    dogname = split(fname, ".")[1]

    # Experimental data
    data = readdlm(fname_exp, ',', Float64, skipstart=1)
    time = data[:,1]
    mid = data[:,5:7]
    tip = data[:,8:10]
    # Simulation data
    data_sim = readdlm(fname_sim, ',', Float64, skipstart=1)
    mid_sim = data_sim[:,2:4]
    tip_sim = data_sim[:,5:7]

    # RMSE
    rmse = ([mid tip] - [mid_sim tip_sim]).^2 |> x -> mean(x, dims=1) |> x -> sqrt.(x) |> mean

    # Plot
    title = plot(title=dogname * " rmse: " * string(round(rmse, digits=3)), grid=false, showaxis=false, ticks=false, bottom_margin=-50Plots.px)
    plt = plot_comparison(mid, tip, mid_sim, tip_sim)
    bigplt = plot(title, plt, layout=@layout([A{0.15h}; B]))
    display(bigplt)
    plts[idx] = bigplt

end
