##### Location of points
function p2(sol, t)
    @inbounds la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb = sol.prob.p
    q1, q2, q3 = sol(t)

    # space123
    # p2x = sol.prob.p.fx(t) + la * cos(q2) * cos(q3)
    # p2y = sol.prob.p.fy(t) + la * sin(q3) * cos(q2)
    # p2z = sol.prob.p.fz(t) - la * sin(q2)

    # body123
    p2x = fx(t) + la * cos(q2) * cos(q3)
    p2y = fy(t) + la * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3))
    p2z = fz(t) + la * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3))

    return p2x, p2y, p2z
end

function p3(sol, t; retp2=false)
    @inbounds la, lb = sol.prob.p
    q1, q2, q3, q4, q5, q6 = sol(t)

    p2x, p2y, p2z = p2(sol, t)

    s1 = sin(q1); s2 = sin(q2); s3 = sin(q3); s4 = sin(q4); s5 = sin(q5); s6 = sin(q6)
    c1 = cos(q1); c2 = cos(q2); c3 = cos(q3); c4 = cos(q4); c5 = cos(q5); c6 = cos(q6)

    # body123
    # p3x = p2x + lb * (cos(q2) * cos(q3) * cos(q5) * cos(q6) - sin(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - sin(q6) * cos(q5) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)))
    # p3y = p2y + lb * (sin(q3) * cos(q2) * cos(q5) * cos(q6) + sin(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) + sin(q6) * cos(q5) * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)))
    # p3z = p2z - lb * (sin(q2) * cos(q5) * cos(q6) + sin(q5) * cos(q1) * cos(q2) - sin(q1) * sin(q6) * cos(q2) * cos(q5))

    # body123
    p3x = p2x + lb * (c2 * c3 * c5 * c6 + s2 * (s4 * s6 - s5 * c4 * c6) - s3 * c2 * (s6 * c4 + s4 * s5 * c6))
    p3y = p2y - lb * (s1 * c2 * (s4 * s6 - s5 * c4 * c6) - c5 * c6 * (s3 * c1 + s1 * s2 * c3) - (s6 * c4 + s4 * s5 * c6) * (c1 * c3 - s1 * s2 * s3))
    p3z = p2z + lb * (c1 * c2 * (s4 * s6 - s5 * c4 * c6) + c5 * c6 * (s1 * s3 - s2 * c1 * c3) + (s1 * c3 + s2 * s3 * c1) * (s6 * c4 + s4 * s5 * c6))
    
    retp2 ? (return p2x, p2y, p2z, p3x, p3y, p3z) : (return p3x, p3y, p3z)
end

p2(sol) = [p2(sol, t) for t ∈ sol.t]
p3(sol; retp2=false) = [p3(sol, t, retp2=retp2) for t ∈ sol.t]

######## Plotting
function plot_model(sol, t; tip=false)
    @inbounds la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb = sol.prob.p
    # px, py, pz = sol.prob.p.f(t)
    px = fx(t)
    py = fy(t)
    pz = fz(t)


    p2x, p2y, p2z, p3x, p3y, p3z = p3(sol, t, retp2=true)

    #  base
    plot([px], [py], [pz], st=:scatter, color=:black, label="")
    # central axes
    plot!([-2,2], [0,0], [0,0], label="", color=:gray, lw=.5)
    plot!([0,0], [-2,2], [0,0], label="", color=:gray, lw=.5)
    plot!([0,0], [0,0], [-2,2], label="", color=:gray, lw=.5)

    # pendulum
    plot!([px, p2x, p3x], [py, p2y, p3y], [pz, p2z, p3z], 
        label="", color=:black, lw=2)# , camera=(0, 30))

    # lines to pendulum tip
    if tip
        plot!([p3x, p3x],[p3y, p3y], [p3z, -2.0], 
              label="", color=:black, ls=:dash)
        plot!([0, p3x], [0, 0], [-1, -1], 
            label="", color=:black, ls=:dash)
        plot!([p3x, p3x], [0, p3y], [-1,-1], 
            label="", color=:black, ls=:dash)
    end

    xlabel!("x"); ylabel!("y"); # zlabel!("z)")
    xlims!(-3, 3)
    ylims!(-3, 3)
    zlims!(-3, 3)
    plot!(aspect_ratio=:equal)

end

function animate_model(sol; fps=24)
    traj = [tup[k] for tup ∈ p3(sol), k ∈ 1:3]
    anim = @animate for t ∈ range(sol.t[begin], sol.t[end], step=1 / fps)
        plot_model(sol, t)
        plot!(traj[:,1], traj[:,2], traj[:,3], label="", ls=:dash, color=:black)
    end

    return gif(anim, "plot.gif", fps=24)
end


##### Energy
function kineticenergy(sol, t)
    la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb = sol.prob.p
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = sol(t)

    # pxp, pyp, pzp = sol.prob.p.fp(t)
    pxp = derivative(fx, t, 1)
    pyp = derivative(fy, t, 1)
    pzp = derivative(fz, t, 1)

    # space 123
    # return 0.5 * ixa * u1^2 + 0.5 * iya * u2^2 + 0.5 * iza * u3^2 + 0.5 * ixb * sin(q5) * u3 * (sin(q5) * u3 - u4 - sin(q6) * cos(q5) * u2 - cos(q5) * cos(q6) * u1) + 0.5 * iyb * u5 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) + 0.5 * izb * u6 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) + 0.5 * iyb * sin(q4) * cos(q5) * u3 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) + 0.5 * izb * cos(q4) * cos(q5) * u3 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) + 0.5 * iyb * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) + 0.5 * izb * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) + 0.5 * mb * (pxp^2 + pyp^2 + pzp^2 + la^2 * u2^2 + la^2 * u3^2 + 2 * la * pzp * sin(q1) * cos(q2) * u3 + 2 * la * pyp * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * u3 + 2 * la * pyp * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) * u2 + lbo^2 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1)^2 + lbo^2 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2)^2 + 2 * la * lbo * cos(q4) * cos(q5) * u2 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) + 2 * la * lbo * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u3 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) + 2 * la * lbo * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u3 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) + 2 * lbo * pzp * (sin(q4) * cos(q1) * cos(q2) * cos(q5) + sin(q2) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) + sin(q1) * cos(q2) * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6))) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) + 2 * lbo * pxp * (sin(q4) * cos(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - cos(q2) * cos(q3) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) - (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3))) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - 2 * la * pzp * cos(q1) * cos(q2) * u2 - 2 * la * pxp * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) * u2 - 2 * la * pxp * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * u3 - 2 * la * lbo * sin(q4) * cos(q5) * u2 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - 2 * lbo * pzp * (cos(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q2) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) - sin(q1) * cos(q2) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) - 2 * lbo * pxp * (cos(q2) * cos(q3) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) + cos(q4) * cos(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) + (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) - 2 * lbo * pyp * (sin(q3) * cos(q2) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) + sin(q4) * cos(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) - (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6))) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - 2 * lbo * pyp * (sin(q3) * cos(q2) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) - cos(q4) * cos(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) - (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1)) - 0.5 * ixb * u4 * (sin(q5) * u3 - u4 - sin(q6) * cos(q5) * u2 - cos(q5) * cos(q6) * u1) - 0.5 * ixb * sin(q6) * cos(q5) * u2 * (sin(q5) * u3 - u4 - sin(q6) * cos(q5) * u2 - cos(q5) * cos(q6) * u1) - 0.5 * ixb * cos(q5) * cos(q6) * u1 * (sin(q5) * u3 - u4 - sin(q6) * cos(q5) * u2 - cos(q5) * cos(q6) * u1) - 0.5 * iyb * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1 * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) - 0.5 * izb * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2 * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - 0.5 * ma * (2 * lao * pzp * cos(q1) * cos(q2) * u2 + 2 * lao * pxp * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) * u2 + 2 * lao * pxp * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * u3 - pxp^2 - pyp^2 - pzp^2 - lao^2 * u2^2 - lao^2 * u3^2 - 2 * lao * pzp * sin(q1) * cos(q2) * u3 - 2 * lao * pyp * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * u3 - 2 * lao * pyp * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) * u2)

    # body123
    return  0.5 * ixa * u1^2 + 0.5 * iya * u2^2 + 0.5 * iza * u3^2 + 0.5 * izb * sin(q4) * cos(q5) * u2 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) + 0.5 * ixb * u4 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) + 0.5 * ixb * cos(q5) * cos(q6) * u1 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) + 0.5 * iyb * sin(q6) * cos(q5) * u1 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + 0.5 * ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) + 0.5 * ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) - 0.5 * izb * u6 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 0.5 * izb * sin(q5) * u1 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 0.5 * izb * cos(q4) * cos(q5) * u3 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 0.5 * iyb * u5 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 0.5 * iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 0.5 * iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 0.5 * ma * (2 * lao * pxp * sin(q2) * u2 + 2 * lao * pxp * sin(q3) * cos(q2) * u3 + 2 * lao * pzp * cos(q1) * cos(q2) * u2 - pxp^2 - pyp^2 - pzp^2 - lao^2 * u2^2 - lao^2 * u3^2 - 2 * lao * pyp * sin(q1) * cos(q2) * u2 - 2 * lao * pzp * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * u3 - 2 * lao * pyp * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * u3) - 0.5 * mb * (2 * la * pxp * sin(q2) * u2 + 2 * la * pxp * sin(q3) * cos(q2) * u3 + 2 * la * pzp * cos(q1) * cos(q2) * u2 + 2 * la * lbo * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u3 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) + 2 * la * lbo * sin(q4) * cos(q5) * u3 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + 2 * la * lbo * cos(q4) * cos(q5) * u2 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + 2 * lbo * pxp * (sin(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - sin(q6) * cos(q2) * cos(q3) * cos(q5) - sin(q3) * cos(q2) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - pxp^2 - pyp^2 - pzp^2 - la^2 * u2^2 - la^2 * u3^2 - 2 * la * pyp * sin(q1) * cos(q2) * u2 - 2 * la * pzp * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * u3 - 2 * la * pyp * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * u3 - lbo^2 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)^2 - 2 * la * lbo * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u2 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - lbo^2 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)^2 - 2 * lbo * pxp * (sin(q2) * cos(q4) * cos(q5) + sin(q5) * cos(q2) * cos(q3) + sin(q3) * sin(q4) * cos(q2) * cos(q5)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 2 * lbo * pyp * (sin(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) + sin(q6) * cos(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 2 * lbo * pzp * (cos(q1) * cos(q2) * cos(q4) * cos(q5) + sin(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - sin(q4) * cos(q5) * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1))) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 2 * lbo * pyp * (sin(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - sin(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q4) * cos(q5) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 2 * lbo * pzp * (sin(q6) * cos(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - cos(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3))

end

function potentialenergy(sol, t)
    la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb = sol.prob.p
    q1, q2, q3, q4, q5, q6 = sol(t)

    pz = evaluate(fz, t)

    # space123
    # return -g * ((ma + mb) * (lao + lbo + pz) - (la * mb + lao * ma) * sin(q2) - lbo * mb * (sin(q2) * cos(q5) * cos(q6) + sin(q5) * cos(q1) * cos(q2) - sin(q1) * sin(q6) * cos(q2) * cos(q5)))
    
    # body123
    return -g * ((ma + mb) * (lao + lbo + pz) + (la * mb + lao * ma) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + lbo * mb * (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))))
    
end

totalenergy(sol, t) = kineticenergy(sol, t) + potentialenergy(sol, t)
kineticenergy(sol) = [kineticenergy(sol, t) for t ∈ sol.t]
potentialenergy(sol) = [potentialenergy(sol, t) for t ∈ sol.t]
totalenergy(sol) = kineticenergy(sol) .+ potentialenergy(sol)

plot_energy(sol) = begin
    plot(potentialenergy(sol), label="pe")
    plot!(kineticenergy(sol), label="ke")
    plot!(totalenergy(sol), label="te")
end

##### Angular momentum about origin
function angmom(sol, t)
    la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb = sol.prob.p
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = sol(t)

    # pxp, pyp, pzp = sol.prob.p.fp(t)
    pxp = derivative(fx, t, 1)
    pyp = derivative(fy, t, 1)
    pzp = derivative(fz, t, 1)

    # space123
    # amomx = lao * ma * (pyp * sin(q2) + pzp * sin(q3) * cos(q2)) + mb * (la * pyp * sin(q2) + la * pzp * sin(q3) * cos(q2) + lbo * pyp * (sin(q2) * cos(q5) * cos(q6) + sin(q5) * cos(q1) * cos(q2) - sin(q1) * sin(q6) * cos(q2) * cos(q5)) + lbo * pzp * (sin(q3) * cos(q2) * cos(q5) * cos(q6) + sin(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) + sin(q6) * cos(q5) * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)))) + cos(q2) * cos(q3) * (ixa * u1 + la * lbo * mb * (sin(q5) * u3 - sin(q6) * cos(q5) * u2)) + (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 + (cos(q2) * cos(q3) * cos(q5) * cos(q6) - sin(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - sin(q6) * cos(q5) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3))) * (ixb * u4 + ixb * sin(q6) * cos(q5) * u2 + ixb * cos(q5) * cos(q6) * u1 - ixb * sin(q5) * u3 - la * lbo * mb * ((sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1))) - (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 - (cos(q2) * cos(q3) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) + cos(q4) * cos(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) + (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (izb * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2 - izb * u6 - izb * cos(q4) * cos(q5) * u3 - izb * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2)) - (sin(q4) * cos(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - cos(q2) * cos(q3) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) - (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3))) * (iyb * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1 - iyb * u5 - iyb * sin(q4) * cos(q5) * u3 - iyb * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1))
    # amomy = sin(q3) * cos(q2) * (ixa * u1 + la * lbo * mb * (sin(q5) * u3 - sin(q6) * cos(q5) * u2)) + (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 + (sin(q3) * cos(q2) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) + sin(q4) * cos(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) - (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6))) * (iyb * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1 - iyb * u5 - iyb * sin(q4) * cos(q5) * u3 - iyb * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1)) + (sin(q3) * cos(q2) * cos(q5) * cos(q6) + sin(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) + sin(q6) * cos(q5) * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3))) * (ixb * u4 + ixb * sin(q6) * cos(q5) * u2 + ixb * cos(q5) * cos(q6) * u1 - ixb * sin(q5) * u3 - la * lbo * mb * ((sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1))) - lao * ma * (pxp * sin(q2) + pzp * cos(q2) * cos(q3)) - mb * (la * pxp * sin(q2) + la * pzp * cos(q2) * cos(q3) + lbo * pxp * (sin(q2) * cos(q5) * cos(q6) + sin(q5) * cos(q1) * cos(q2) - sin(q1) * sin(q6) * cos(q2) * cos(q5)) + lbo * pzp * (cos(q2) * cos(q3) * cos(q5) * cos(q6) - sin(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - sin(q6) * cos(q5) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)))) - (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 - (sin(q3) * cos(q2) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) - cos(q4) * cos(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) - (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (izb * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2 - izb * u6 - izb * cos(q4) * cos(q5) * u3 - izb * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2))
    # amomz = sin(q1) * cos(q2) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 + cos(q1) * cos(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 - lao * ma * cos(q2) * (pxp * sin(q3) - pyp * cos(q3)) - mb * (la * pxp * sin(q3) * cos(q2) + lbo * pxp * (sin(q3) * cos(q2) * cos(q5) * cos(q6) + sin(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) + sin(q6) * cos(q5) * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3))) - la * pyp * cos(q2) * cos(q3) - lbo * pyp * (cos(q2) * cos(q3) * cos(q5) * cos(q6) - sin(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - sin(q6) * cos(q5) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)))) - sin(q2) * (ixa * u1 + la * lbo * mb * (sin(q5) * u3 - sin(q6) * cos(q5) * u2)) - (sin(q4) * cos(q1) * cos(q2) * cos(q5) + sin(q2) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) + sin(q1) * cos(q2) * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6))) * (iyb * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1 - iyb * u5 - iyb * sin(q4) * cos(q5) * u3 - iyb * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1)) - (cos(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q2) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) - sin(q1) * cos(q2) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (izb * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2 - izb * u6 - izb * cos(q4) * cos(q5) * u3 - izb * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2)) - (sin(q2) * cos(q5) * cos(q6) + sin(q5) * cos(q1) * cos(q2) - sin(q1) * sin(q6) * cos(q2) * cos(q5)) * (ixb * u4 + ixb * sin(q6) * cos(q5) * u2 + ixb * cos(q5) * cos(q6) * u1 - ixb * sin(q5) * u3 - la * lbo * mb * ((sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1)))

    # body123
    amomx = lao * ma * (pzp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - pyp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3))) + mb * (la * pzp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - la * pyp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - lbo * pyp * (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - lbo * pzp * (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)))) + sin(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 + cos(q2) * cos(q3) * (ixa * u1 - la * lbo * mb * ((sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3)) + (sin(q2) * cos(q4) * cos(q5) + sin(q5) * cos(q2) * cos(q3) + sin(q3) * sin(q4) * cos(q2) * cos(q5)) * (izb * u6 + izb * sin(q5) * u1 + izb * cos(q4) * cos(q5) * u3 - izb * sin(q4) * cos(q5) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)) + (sin(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - sin(q6) * cos(q2) * cos(q3) * cos(q5) - sin(q3) * cos(q2) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (iyb * u5 + iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 + iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 - iyb * sin(q6) * cos(q5) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)) + (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) * (ixb * u4 + ixb * cos(q5) * cos(q6) * u1 + ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 + la * lbo * mb * (sin(q5) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - sin(q6) * cos(q5) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2))) - sin(q3) * cos(q2) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2
    amomy = (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 + (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) * (ixa * u1 - la * lbo * mb * ((sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3)) + (sin(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - sin(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q4) * cos(q5) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (izb * u6 + izb * sin(q5) * u1 + izb * cos(q4) * cos(q5) * u3 - izb * sin(q4) * cos(q5) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)) - lao * ma * (pzp * cos(q2) * cos(q3) - pxp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3))) - mb * (la * pzp * cos(q2) * cos(q3) + lbo * pzp * (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - la * pxp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - lbo * pxp * (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)))) - sin(q1) * cos(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 - (sin(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) + sin(q6) * cos(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (iyb * u5 + iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 + iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 - iyb * sin(q6) * cos(q5) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)) - (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (ixb * u4 + ixb * cos(q5) * cos(q6) * u1 + ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 + la * lbo * mb * (sin(q5) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - sin(q6) * cos(q5) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)))
    amomz = lao * ma * (pyp * cos(q2) * cos(q3) - pxp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3))) + cos(q1) * cos(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 + (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) * (ixa * u1 - la * lbo * mb * ((sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3)) + (cos(q1) * cos(q2) * cos(q4) * cos(q5) + sin(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - sin(q4) * cos(q5) * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1))) * (izb * u6 + izb * sin(q5) * u1 + izb * cos(q4) * cos(q5) * u3 - izb * sin(q4) * cos(q5) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)) + (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) * (ixb * u4 + ixb * cos(q5) * cos(q6) * u1 + ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 + la * lbo * mb * (sin(q5) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - sin(q6) * cos(q5) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2))) - mb * (la * pxp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - la * pyp * cos(q2) * cos(q3) - lbo * pyp * (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - lbo * pxp * (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)))) - (sin(q6) * cos(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - cos(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (iyb * u5 + iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 + iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 - iyb * sin(q6) * cos(q5) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2))

    return amomx, amomy, amomz
end

angmom(sol) = [tup[k] for tup ∈ [angmom(sol, t) for t ∈ sol.t], k ∈ 1:3]
