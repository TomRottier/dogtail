using Dierckx, MAT
# Extract data from .mat files and fit quintic spline to base position
# Returns splines for each direction as function of time and position data for base, midpoint and tip of tail
function getdata(fname)
    # Load data
    # fname = "MPIData/MPI-Advanced-Kay0040.mat";
    data = matread(fname) |> values |> collect; data = data[1]

    # Get tail marker data
    markers = data["Trajectories"]["Labeled"]
    mdata = markers["Data"] |> x -> permutedims(x, [3,2,1])
    base_marker = "SchwAns"; mid_marker = "Schw"; tip_marker = "SchwSpi"
    mnames = markers["Labels"] |> vec

    base = mdata[:, 1:3, base_marker .== mnames] ./ 1000 # Get in m
    mid  = mdata[:, 1:3, mid_marker .== mnames]  ./ 1000 
    tip  = mdata[:, 1:3, tip_marker .== mnames]  ./ 1000 

    # Plot to check
    # anim = @animate for t ∈ 1:size(base, 1)
    #     plot([base[t,1],mid[t,1],tip[t,1]], 
    #     [base[t,2],mid[t,2],tip[t,2]], 
    #     [base[t,3],mid[t,3],tip[t,3]], lw=2, color=:black, label="" )
    #     plot!([base[t,1]], [base[t,2]], [base[t,3]], markershape=:circle, markersize=4, color=:black, label="")
    #     mins = minimum([base; mid; tip], dims=1); maxs = maximum([base; mid; tip], dims=1)
    #     xlims!(mins[1], maxs[1])
    #     ylims!(mins[2], maxs[2])
    #     zlims!(mins[3], maxs[3])
    #     plot!(aspect_ratio=:equal)
    # end

    # gif(anim, "recorded.gif")


    # Fit quintic spline to base data
    hz = data["FrameRate"]; n = data["Frames"]
    time = range(0, step=1 / hz, length=Int(n))
    splx = Spline1D(time, base[:,1]) 
    sply = Spline1D(time, base[:,2]) 
    splz = Spline1D(time, base[:,3])

    return splx, sply, splz, base, mid, tip
end