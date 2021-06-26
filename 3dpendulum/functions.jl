getx(sol) = @. sol.prob.p[6] * cos(sol[2,:]) * cos(sol[3,:])
gety(sol) = @. sol.prob.p[6] * sin(sol[3,:]) * cos(sol[2,:])
getz(sol) = @. -sol.prob.p[6] * sin(sol[2,:])
getx(sol,t) = sol.prob.p[6] * cos(sol(t)[2]) * cos(sol(t)[3])
gety(sol,t) = sol.prob.p[6] * sin(sol(t)[3]) * cos(sol(t)[2])
getz(sol,t) = -sol.prob.p[6] * sin(sol(t)[2])

function plot_pendulum(sol, t; tip=false)
    plot([0], [0], [0], st=:scatter, color=:black, label="")
    plot!([0,0], [-2,2], [0,0], label="", color=:gray, lw=.5)
    plot!([-2,2], [0,0], [0,0], label="", color=:gray, lw=.5)
    plot!([0,0], [0,0], [-2,2], label="", color=:gray, lw=.5)
    plot!([0, getx(sol, t)], 
         [0, gety(sol, t)], 
         [0, getz(sol, t)],
         label="", color=:black, lw=2)# , camera=(0, 30))
    if tip
        plot!([getx(sol, t),getx(sol, t)],
              [gety(sol, t),gety(sol, t)],
              [getz(sol, t), -2.0], 
              label="", color=:black, ls=:dash)
        plot!([0, getx(sol, t)], [0,0], [-1,-1], 
            label="", color=:black, ls=:dash)
        plot!([getx(sol, t), getx(sol, t)], [0,gety(sol, t)], [-1,-1], 
            label="", color=:black, ls=:dash)
    end
    # println(getx(sol,t), gety(sol,t), getz(sol,t))

    xlabel!("x"); ylabel!("y"); # zlabel!("z)")
    xlims!(-1, 1)
    ylims!(-1, 1)
    zlims!(-1, 1)
    plot!(aspect_ratio=:equal)

end
function animate_pendulum(sol;fps=24)
    px = getx(sol); py = gety(sol); pz = getz(sol)
    anim = @animate for t ∈ range(sol.t[begin], sol.t[end], step=1 / fps)
        plot_pendulum(sol, t)
        plot!(px, py, pz, label="", ls=:dash, color=:black)
    end

    return gif(anim, "plot.gif", fps=24)
end


###### Energy
# ke = 0.5d0*ix*u4^2 + 0.5d0*iy*u5^2 + 0.5d0*iz*u6^2 + 0.5d0*m*lao^2*(u5^2+u6^2)
# pe = -g*m*(q3-lao*sin(q5))
# te = pe + ke
function kineticenergy(sol, t)
    ma, ix, iy, iz, lao, la, g = sol.prob.p
    q4, q5, q6, u4, u5, u6 = sol(t)

    return 0.5 * ix * u4^2 + 0.5 * iy * u5^2 + 0.5 * iz * u6^2 + 0.5 * m * lao^2 * (u5^2 + u6^2)
end

function potentialenergy(sol, t)
    ma, ix, iy, iz, lao, la, g = sol.prob.p
    q4, q5, q6, u4, u5, u6 = sol(t)

    return -g * m * (- lao * sin(q5))
end

totalenergy(sol, t) = kineticenergy(sol, t) + potentialenergy(sol, t)
kineticenergy(sol) = [kineticenergy(sol, t) for t ∈ sol.t]
potentialenergy(sol) = [potentialenergy(sol, t) for t ∈ sol.t]
totalenergy(sol) = kineticenergy(sol) .+ potentialenergy(sol)
