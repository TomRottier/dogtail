# Cost function to optimise
function cost(pin)
    # Remake problem with new parameters
    pvariable = (ka = pin[1], kb = pin[2], ba = pin[3], bb = pin[4], ma = pin[5], mb = pin[6])
    pnew = merge(p, pvariable)
    # newprob = remake(prob, p=pnew)
    newprob = ODEProblem(eom!, u₀, tspan, pnew)
    

    # Solve
    sol = solve(newprob, Tsit5(), abstol=1e-8, reltol=1e-7)

    # Compare simulation output to measured data
    mid_sim = [tup[j] for tup ∈ [p2(sol, t) for t ∈ times], j ∈ 1:3]
    tip_sim = [tup[j] for tup ∈ [p3(sol, t) for t ∈ times], j ∈ 1:3]

    # Average RMSE for each direction
    mid_rmse = sqrt.((mean((mid - mid_sim).^2, dims=1))) |> mean
    tip_rmse = sqrt.((mean((tip - tip_sim).^2, dims=1))) |> mean

    return mid_rmse^2 + tip_rmse^2
end
