# Extract data from .mat files and fit quintic spline to base position
# Returns splines for each direction as function of time, position data for base, midpoint and tip of tail, orientations of tail segments and average length of tails
using Dierckx, MAT, DSP

# Filter file with specified cutoff and save output to csv
function getdata(fname, fc; file=false, showplot=true)
    # Load data
    # fname = "MPIData/MPI-Advanced-Kay0040.mat";
    data = matread(fname) |> values |> collect; data = data[1]
    fs = data["FrameRate"]

    # Get tail marker data
    markers = data["Trajectories"]["Labeled"]
    mdata = markers["Data"] |> x -> permutedims(x, [3,2,1]) 
    base_marker = "SchwAns"; mid_marker = "Schw"; tip_marker = "SchwSpi"
    rHip_marker = "RISG"; lHip_marker = "LISG"
    mnames = markers["Labels"] |> vec

    # Get individual marker data in m
    base_raw = mdata[:, 1:3, base_marker .== mnames][:,:,1] ./ 1000
    mid_raw  = mdata[:, 1:3, mid_marker .== mnames][:,:,1]  ./ 1000
    tip_raw  = mdata[:, 1:3, tip_marker .== mnames][:,:,1]  ./ 1000
    rHip_raw = mdata[:, 1:3, rHip_marker .== mnames][:,:,1] ./ 1000
    lHip_raw = mdata[:, 1:3, lHip_marker .== mnames][:,:,1] ./ 1000

    # Fit interpolating spline to data then evaluate on same time base to remove NaN's
    n = size(mdata, 1)
    time = range(0, step=1 / fs, length=n) |> collect   # Common time base
    time_end = time[end]
    dirty = [base_raw, mid_raw, tip_raw, rHip_raw, lHip_raw]        # Collect data potentially with NaN's
    clean = Vector(undef, 5)
    for i ∈ 1:5
        _data = dirty[i] |> x -> x[ .!isnan.(x) ] |> x -> reshape(x, :, 3)
        _n = size(_data, 1)
        time_dirty = range(0, time_end, length=_n) |> collect
        splx = Spline1D(time_dirty, _data[:,1], k=5)
        sply = Spline1D(time_dirty, _data[:,2], k=5)
        splz = Spline1D(time_dirty, _data[:,3], k=5)

        data_cleaned = [evaluate(splx, time) evaluate(sply, time) evaluate(splz, time)]
        clean[i] = data_cleaned
    end

    base, mid, tip, rHip, lHip = clean      # Data with NaN's removed, on common time base
    
    # Filter data; low pass 2nd order Butterworth filter (passed twice so effectively 4th order, fc not adjusted for double pass so actual fc lower than fc)
    f = digitalfilter(Lowpass(fc, fs=fs), Butterworth(2))
    base_f = filtfilt(f, base); mid_f = filtfilt(f, mid); tip_f = filtfilt(f, tip)
    rHip_f = filtfilt(f, rHip); lHip_f = filtfilt(f, lHip); 
    hip_f = 0.5 * (rHip_f + lHip_f)

    # Calculate orientation of hip frame
    f_ang = getorientation(base_f - hip_f, [[1.,0.,0.] for _ ∈ 1:size(base_f, 1)])
    splx = Spline1D(time, f_ang[:,1], k=5)
    sply = Spline1D(time, f_ang[:,2], k=5)
    splz = Spline1D(time, f_ang[:,3], k=5)

    # Plot output
    if showplot
        # Plot data
        # plt = plot(base_raw, label=["x_raw" "y_raw" "z_raw"], title="base")
        # plot!(base, label=["x" "y" "z"])
        # plot!(base_f, label=["x_filt" "y_filt" "z_filt"])
        # display(plt)
        # plt = plot(mid_raw, label=["x_raw" "y_raw" "z_raw"], title="mid")
        # plot!(mid, label=["x" "y" "z"])
        # plot!(mid_f, label=["x_filt" "y_filt" "z_filt"])
        # display(plt)
        # plt = plot(tip_raw, label=["x_raw" "y_raw" "z_raw"], title="tip")
        # plot!(tip, label=["x" "y" "z"])
        # plot!(tip_f, label=["x_filt" "y_filt" "z_filt"])
        # display(plt)

        # Plot spline 2nd derivative
        plt = plot(derivative(splx, time, 2), label="x_acc", title=fname)
        plot!(derivative(sply, time, 2), label="y_acc")
        plot!(derivative(splz, time, 2), label="z_acc")
        display(plt)
    end

    if file
        name = split(fname, "-") |> x -> split(x[3], ".mat")[1]
        out = [time f_ang base_f mid_f tip_f]
        header = ["time" "hipX" "hipY" "hipZ" "baseX" "baseY" "baseZ" "midX" "midY" "midZ" "tipX" "tipY" "tipZ"]
        folder = "Data//"
        open(folder * name * ".csv", "w") do io
            writedlm(io, [header; out], ',')
        end
    end
end

# Read data from csv
function getdata(fname)
    # Load data
    data = readdlm(fname, ',', skipstart=1)
    time = data[:,1]
    f_ang = data[:,2:4]
    base = data[:,5:7]
    mid = data[:,8:10]
    tip = data[:,11:13]
    
    # Fit quintic interpolating splines
    splx = Spline1D(time, base[:,1], k=5) # Spline1D(time, repeat([0.0], size(time, 1)), k=5)  
    sply = Spline1D(time, base[:,2], k=5) # Spline1D(time, repeat([0.0], size(time, 1)), k=5)  
    splz = Spline1D(time, base[:,3], k=5) # Spline1D(time, repeat([0.0], size(time, 1)), k=5)  
    splθ₁ = Spline1D(time, f_ang[:,1], k=5)  # Spline1D(time, repeat([0.0], size(time, 1)), k=5) 
    splθ₂ = Spline1D(time, f_ang[:,2], k=5)  # Spline1D(time, repeat([0.0], size(time, 1)), k=5) 
    splθ₃ = Spline1D(time, f_ang[:,3], k=5)  # Spline1D(time, repeat([0.0], size(time, 1)), k=5) 

    # Orientation angles body123 of body a relative to frame f (hip frame)
    ref1 = [RotXYZ(vec...) * [1,0,0] |> Array for vec ∈ eachrow(f_ang)]
    ref2 = [mid[i,:] - base[i,:] for i ∈ 1:size(base, 1)]
    or1 = getorientation(mid - base, ref1)
    or2 = getorientation(tip - mid, ref2)

    # Extract Euler angles
    q1 = or1[:,1]; q2 = or1[:,2]; q3 = or1[:,3]
    q4 = or2[:,1]; q5 = or2[:,2]; q6 = or2[:,3]

    # Fit splines to Euler angles
    spl_q1 = Spline1D(time, q1, k=5); spl_q2 = Spline1D(time, q2, k=5); spl_q3 = Spline1D(time, q3, k=5)
    spl_q4 = Spline1D(time, q4, k=5); spl_q5 = Spline1D(time, q5, k=5); spl_q6 = Spline1D(time, q6, k=5)

    # First derivative of Euler angles
    q1p = derivative(spl_q1, time, 1); q2p = derivative(spl_q2, time, 1); q3p = derivative(spl_q3, time, 1)
    q4p = derivative(spl_q4, time, 1); q5p = derivative(spl_q5, time, 1); q6p = derivative(spl_q6, time, 1)
    
    # Convert to generalised speeds
    u1 = @. sin(q3) * q2p + cos(q2) * cos(q3) * q1p
    u2 = @. cos(q3) * q2p - sin(q3) * cos(q2) * q1p
    u3 = @. q3p + q1p * sin(q2)
    u4 = @. sin(q6) * q5p + cos(q5) * cos(q6) * q4p
    u5 = @. cos(q6) * q5p - sin(q6) * cos(q5) * q4p
    u6 = @. q6p + q4p * sin(q5)

    initconds = [q1[1],q2[1],q3[1],q4[1],q5[1],q6[1],u1[1],u2[1],u3[1],u4[1],u5[1],u6[1]]

    # Segment lengths; removing NaN's
    la = sqrt.(sum((mid - base).^2, dims=2)) |> mean
    lb = sqrt.(sum((tip - mid).^2, dims=2))  |> mean

    # Time span
    tspan = (0, time[end])

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


    return splx, sply, splz, splθ₁, splθ₂, splθ₃, time, base, mid, tip, initconds, la, lb, tspan
end

# Process each file individually 
# fnames = readdir("MPIData/")
# for fname in fnames
#     getdata("MPIData/" *  fname, 8, file=false, showplot=true)
# end