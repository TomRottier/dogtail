# Extract data from .mat files and fit quintic spline to base position
# Returns splines for each direction as function of time, position data for base, midpoint and tip of tail, orientations of tail segments and average length of tails
using Dierckx, MAT, DSP, DelimitedFiles

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

    # gif(anim, "recorded.gif", fps=48)


    # Fit quintic smoothing spline to data then evaluate across common time to remove and NaN's
    hz = data["FrameRate"]; n = size(mdata, 1)
    time = range(0, step=1 / hz, length=n) |> collect
    time_end = time[end]
    dirty = [base, mid, tip]
    clean = Vector(undef, 3)
    for i ∈ 1:3
        _data = dirty[i] |> x -> x[ .!isnan.(x) ] |> x -> reshape(x, :, 3)
        n = size(_data, 1)
        time_dirty = range(0, time_end, length=n) |> collect
        splx = Spline1D(time_dirty, _data[:,1], k=5, s=0.0001)
        sply = Spline1D(time_dirty, _data[:,2], k=5, s=0.0001)
        splz = Spline1D(time_dirty, _data[:,3], k=5, s=0.0001)

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

# Filter file with specified cutoff and save output to csv
function getdata(fname, fc; file=false, plot=true)
    # Load data
    # fname = "MPIData/MPI-Advanced-Kay0040.mat";
    data = matread(fname) |> values |> collect; data = data[1]
    fs = data["FrameRate"]

    # Get tail marker data
    markers = data["Trajectories"]["Labeled"]
    mdata = markers["Data"] |> x -> permutedims(x, [3,2,1]) 
    base_marker = "SchwAns"; mid_marker = "Schw"; tip_marker = "SchwSpi"
    mnames = markers["Labels"] |> vec

    # Get individual marker data in m
    base_raw = mdata[:, 1:3, base_marker .== mnames][:,:,1] ./ 1000
    mid_raw  = mdata[:, 1:3, mid_marker .== mnames][:,:,1]  ./ 1000
    tip_raw  = mdata[:, 1:3, tip_marker .== mnames][:,:,1]  ./ 1000

    # Fit interpolating spline to data then evaluate on same time base to remove NaN's
    n = size(mdata, 1)
    time = range(0, step=1 / hz, length=n) |> collect   # Common time base
    time_end = time[end]
    dirty = [base_raw, mid_raw, tip_raw]        # Collect data potentially with NaN's
    clean = Vector(undef, 3)
    for i ∈ 1:3
        _data = dirty[i] |> x -> x[ .!isnan.(x) ] |> x -> reshape(x, :, 3)
        _n = size(_data, 1)
        time_dirty = range(0, time_end, length=_n) |> collect
        splx = Spline1D(time_dirty, _data[:,1], k=5)
        sply = Spline1D(time_dirty, _data[:,2], k=5)
        splz = Spline1D(time_dirty, _data[:,3], k=5)

        data_cleaned = [evaluate(splx, time) evaluate(sply, time) evaluate(splz, time)]
        clean[i] = data_cleaned
    end

    base, mid, tip = clean      # Data with NaN's removed, on common time base
    
    # Filter data; low pass 2nd order Butterworth filter (passed twice so effectively 4th order, fc not adjusted for double pass so actual fc lower than fc)
    f = digitalfilter(Lowpass(fc, fs=fs), Butterworth(2))
    base_f = filtfilt(f, base); mid_f = filtfilt(f, mid); tip_f = filtfilt(f, tip)

    # Fit quintic interpolating spline to filtered base data
    splx = Spline1D(time, base_f[:,1], k=5)
    sply = Spline1D(time, base_f[:,2], k=5)
    splz = Spline1D(time, base_f[:,3], k=5)

    # Plot output
    if plot
        # Plot data
        plt = plot(base_raw, label=["x_raw" "y_raw" "z_raw"], title="base")
        plot!(base, label=["x" "y" "z"])
        plot!(base_f, label=["x_filt" "y_filt" "z_filt"])
        display(plt)
        plt = plot(mid_raw, label=["x_raw" "y_raw" "z_raw"], title="mid")
        plot!(mid, label=["x" "y" "z"])
        plot!(mid_f, label=["x_filt" "y_filt" "z_filt"])
        display(plt)
        plt = plot(tip_raw, label=["x_raw" "y_raw" "z_raw"], title="tip")
        plot!(tip, label=["x" "y" "z"])
        plot!(tip_f, label=["x_filt" "y_filt" "z_filt"])
        display(plt)

        # Plot spline 2nd derivative
        plt = plot(derivative(splx, time, 2), label="x_acc")
        plot!(derivative(sply, time, 2), label="y_acc")
        plot!(derivative(splz, time, 2), label="z_acc")
        display(plt)
    end

    if file
        name = split(fname, "-") |> x -> split(x[3], ".mat")[1]
        out = [time base mid tip]
        header = ["time" "baseX" "baseY" "baseZ" "midX" "midY" "midZ" "tipX" "tipY" "tipZ"]
        folder = "Data//"
        open(folder * name * ".csv", "w") do io
            writedlm(io, [header; out], ',')
        end
        
    end



    
end

