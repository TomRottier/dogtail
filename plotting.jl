function animate_comparison(base, mid, tip, mid_sim, tip_sim)
    base_sim = copy(base)
    mins = minimum([base; mid; tip; base_sim; mid_sim; tip_sim], dims=1); maxs = maximum([base; mid; tip; base_sim; mid_sim; tip_sim], dims=1)

    anim = @animate for (index, t) ∈ enumerate(times)
        # Trajectory
        plot(base[:,1], base[:,2], base[:,3], label="", ls=:dash, color=:black)
        # Base
        plot!([base[index,1]], [base[index,2]], [base[index,3]], label="", st=:scatter, color=:blue)
        # Actual tail
        plot!([base[index,1], mid[index,1], tip[index,1]],
             [base[index,2], mid[index,2], tip[index,2]],
             [base[index,3], mid[index,3], tip[index,3]],
             color=:black, lw=:2, label="")
        # Simulated tail
        plot!([base_sim[index,1], mid_sim[index,1], tip_sim[index,1]],
             [base_sim[index,2], mid_sim[index,2], tip_sim[index,2]],
             [base_sim[index,3], mid_sim[index,3], tip_sim[index,3]],
             color=:red, lw=:2, label="")
        # Axis limits
        xlims!(mins[1], maxs[1])
        ylims!(mins[2], maxs[2])
        zlims!(mins[3], maxs[3])
        plot!(aspect_ratio=:equal)
        # Camera angle
        # plot!(camera=(90,60))
    end
    
    return gif(anim, "comparison.gif", fps=24)
end

function plot_comparison(mid, tip, mid_sim, tip_sim)

    plt1 = plot(mid, label=["x" "y" "z"]); plot!(mid_sim, label=["x_sim" "y_sim" "z_sim"]); plot!(title="mid", legend=:none);
    plt2 = plot(tip, label=["x" "y" "z"]); plot!(tip_sim, label=["x_sim" "y_sim" "z_sim"]); plot!(title="tip");
    return plot(plt1, plt2, layout=(1, 2), size=(1000, 400))

end