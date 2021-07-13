# Cost function to optimise
function cost(pin, p, prob, times, mid, tip)
    # Remake problem with new parameters
    update_parameters!(pin, p)
    # newprob = ODEProblem{true}(eom!, u₀, tspan, p)
    newprob = remake(prob; p=p)

    # Solve
    sol = solve(newprob, Tsit5(), abstol=1e-7, reltol=1e-7, saveat=times)
    sol.retcode == :MaxIters && return 1.0

    # Compare simulation output to measured data
    mid_sim::Matrix{Float64} = [tup[i] for tup ∈ p2(sol), i ∈ 1:3]
    tip_sim::Matrix{Float64} = [tup[i] for tup ∈ p3(sol), i ∈ 1:3]

    # Average RMSE for each direction
    mid_rmse = sqrt.((mean((mid::Matrix{Float64} - mid_sim).^2, dims=1))) |> mean
    tip_rmse = sqrt.((mean((tip::Matrix{Float64} - tip_sim).^2, dims=1))) |> mean

    return (mid_rmse + tip_rmse) * 0.5
end

# Returns simulated values
function cost(pin, p, u₀, tspan, times)
    # Remake problem with new parameters
    update_parameters!(pin, p)
    newprob = ODEProblem{true}(eom!, u₀, tspan, p)
    
    # Solve
    sol = solve(newprob, Tsit5(), abstol=1e-7, reltol=1e-7, saveat=times)
    sol.retcode == :MaxIters && return 1.0

    # Compare simulation output to measured data
    mid_sim::Matrix{Float64} = [tup[i] for tup ∈ p2(sol), i ∈ 1:3]
    tip_sim::Matrix{Float64} = [tup[i] for tup ∈ p3(sol), i ∈ 1:3]

    return mid_sim, tip_sim
end
