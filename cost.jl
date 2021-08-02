# Cost function to optimise
function cost(pin, p, prob, times, mid, tip)
    # Remake problem with new parameters
    update_parameters!(p, pin)
    # newprob = ODEProblem{true}(eom!, u₀, tspan, p)
    newprob = remake(prob; p=p)

    # Solve
    sol = solve(newprob, Tsit5(), abstol=1e-5, reltol=1e-5, saveat=times)
    sol.retcode == :MaxIters  && (println(pin),return 1.0)
    sol.retcode == :DtLessThanMin && return 2.0

    # Compare simulation output to measured data
    mid_sim::Matrix{Float64} = [p2x(sol) p2y(sol) p2z(sol)]
    tip_sim::Matrix{Float64} = [p3x(sol) p3y(sol) p3z(sol)]
    
    # Average RMSE for each direction
    mid_rmse = sqrt.((mean((mid::Matrix{Float64} - mid_sim).^2, dims=1))) |> mean
    tip_rmse = sqrt.((mean((tip::Matrix{Float64} - tip_sim).^2, dims=1))) |> mean

    return 0.8 * mid_rmse + 0.2 * tip_rmse
end

# Returns simulated values
function cost(pin, p, prob, times)
    # Remake problem with new parameters
    update_parameters!(p, pin)
    newprob = remake(prob; p=p)
    
    # Solve
    sol = solve(newprob, Tsit5(), abstol=1e-7, reltol=1e-7, saveat=times)
    sol.retcode == :MaxIters && return 1.0

    # Compare simulation output to measured data
    mid_sim::Matrix{Float64} = [p2x(sol) p2y(sol) p2z(sol)]
    tip_sim::Matrix{Float64} = [p3x(sol) p3y(sol) p3z(sol)]

    return mid_sim, tip_sim
end
