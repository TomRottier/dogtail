# Extract data from .mat files and fit quintic spline to base position
# Returns splines for each direction as function of time, position data for base, midpoint and tip of tail, orientations of tail segments and average length of tails
using Dierckx, MAT

function getdata(fname)
    # Load data
    # fname = "MPIData/MPI-Advanced-Kay0040.mat";
    data = matread(fname) |> values |> collect; data = data[1]

    # Get tail marker data
    markers = data["Trajectories"]["Labeled"]
    mdata = markers["Data"] |> x -> permutedims(x, [3,2,1]) 
    base_marker = "SchwAns"; mid_marker = "Schw"; tip_marker = "SchwSpi"
    mnames = markers["Labels"] |> vec

    # Ge individual marker data in m
    base = mdata[:, 1:3, base_marker .== mnames][:,:,1] ./ 1000
    mid  = mdata[:, 1:3, mid_marker .== mnames][:,:,1]  ./ 1000
    tip  = mdata[:, 1:3, tip_marker .== mnames][:,:,1]  ./ 1000

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


    # Fit quintic spline to data then evaluate across common time to remove and NaN's
    hz = data["FrameRate"]; n = size(mdata, 1)
    time = range(0, step=1 / hz, length=n)
    dirty = [base, mid, tip]
    clean = Vector(undef, 3)
    for i ∈ 1:3
        data = dirty[i] |> x -> x[ .!isnan.(x) ] |> x -> reshape(x, :, 3)
        n = size(data, 1)
        time_dirty = range(0, step=1/hz, length=n)
        splx = Spline1D(time_dirty, data[:,1], k=5)
        sply = Spline1D(time_dirty, data[:,2], k=5)
        splz = Spline1D(time_dirty, data[:,3], k=5)

        data_cleaned = [evaluate(splx, time) evaluate(sply, time) evaluate(splz, time)]
        clean[i] = data_cleaned
    end

    base, mid, tip = clean

    # Output splines
    splx = Spline1D(time, base[:,1], k=5)
    sply = Spline1D(time, base[:,2], k=5)
    splz = Spline1D(time, base[:,3], k=5)


    # Orientation angles body123
    orientations = getorientation(time, base, mid, tip)

    # Segment lengths; removing NaN's
    la = sqrt.(sum((mid - base).^2, dims=2)) |> mean
    lb = sqrt.(sum((tip - mid).^2, dims=2))  |> mean

    # Time span
    tspan = (0, time[end])

    return splx, sply, splz, base, mid, tip, orientations, la, lb, tspan
end