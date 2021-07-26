##### Location of points
function p2(sol, t)
    @unpack la, fx, fy, fz = sol.prob.p
    q1, q2, q3 = sol(t)

    px = fx(t); py = fy(t); pz = fz(t)

    s1 = sin(q1); s2 = sin(q2); s3 = sin(q3);
    c1 = cos(q1); c2 = cos(q2); c3 = cos(q3);
    
    p2x = px + la * c2 * c3
    p2y = py + la * (s3 * c1 + s1 * s2 * c3)
    p2z = pz + la * (s1 * s3 - s2 * c1 * c3)
        
    return p2x, p2y, p2z
end

function p3(sol, t; retp2=false)
    @unpack la, lb, fx, fy, fz = sol.prob.p
    q1, q2, q3, q4, q5, q6 = sol(t)

    px = fx(t); py = fy(t); pz = fz(t)

    retp2 && (p2x, p2y, p2z = p2(sol, t))

    s1 = sin(q1); s2 = sin(q2); s3 = sin(q3); s4 = sin(q4); s5 = sin(q5); s6 = sin(q6)
    c1 = cos(q1); c2 = cos(q2); c3 = cos(q3); c4 = cos(q4); c5 = cos(q5); c6 = cos(q6)

    p3x = px + la * c2 * c3 + lb * (c2 * c3 * c5 * c6 + s2 * (s4 * s6 - s5 * c4 * c6) - s3 * c2 * (s6 * c4 + s4 * s5 * c6))
    p3y = py + la * (s3 * c1 + s1 * s2 * c3) - lb * (s1 * c2 * (s4 * s6 - s5 * c4 * c6) - c5 * c6 * (s3 * c1 + s1 * s2 * c3) - (s6 * c4 + s4 * s5 * c6) * (c1 * c3 - s1 * s2 * s3))
    p3z = pz + la * (s1 * s3 - s2 * c1 * c3) + lb * (c1 * c2 * (s4 * s6 - s5 * c4 * c6) + c5 * c6 * (s1 * s3 - s2 * c1 * c3) + (s1 * c3 + s2 * s3 * c1) * (s6 * c4 + s4 * s5 * c6))
            
    retp2 ? (return p2x, p2y, p2z, p3x, p3y, p3z) : (return p3x, p3y, p3z)
end

p2(sol) = [p2(sol, t) for t ∈ sol.t]
p3(sol; retp2=false) = [p3(sol, t, retp2=retp2) for t ∈ sol.t]


##### Energy
function kineticenergy(sol, t)
    @unpack la, lao, lbo, ixa, ixb, iya, iyb, iza, izb, fx, fy, fz, ma, mb = sol.prob.p
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = sol(t)

    # pxp, pyp, pzp = sol.prob.p.fp(t)
    pxp = derivative(fx, t, 1)
    pyp = derivative(fy, t, 1)
    pzp = derivative(fz, t, 1)

    return 0.5 * ixa * u1^2 + 0.5 * iya * u2^2 + 0.5 * iza * u3^2 + 0.5 * izb * sin(q4) * cos(q5) * u2 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) + 0.5 * ixb * u4 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) + 0.5 * ixb * cos(q5) * cos(q6) * u1 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) + 0.5 * iyb * sin(q6) * cos(q5) * u1 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + 0.5 * ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) + 0.5 * ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 * (u4 + cos(q5) * cos(q6) * u1 + (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3) - 0.5 * izb * u6 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 0.5 * izb * sin(q5) * u1 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 0.5 * izb * cos(q4) * cos(q5) * u3 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 0.5 * iyb * u5 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 0.5 * iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 0.5 * iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 0.5 * ma * (2 * lao * pxp * sin(q2) * u2 + 2 * lao * pxp * sin(q3) * cos(q2) * u3 + 2 * lao * pzp * cos(q1) * cos(q2) * u2 - pxp^2 - pyp^2 - pzp^2 - lao^2 * u2^2 - lao^2 * u3^2 - 2 * lao * pyp * sin(q1) * cos(q2) * u2 - 2 * lao * pzp * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * u3 - 2 * lao * pyp * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * u3) - 0.5 * mb * (2 * la * pxp * sin(q2) * u2 + 2 * la * pxp * sin(q3) * cos(q2) * u3 + 2 * la * pzp * cos(q1) * cos(q2) * u2 + 2 * la * lbo * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u3 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) + 2 * la * lbo * sin(q4) * cos(q5) * u3 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + 2 * la * lbo * cos(q4) * cos(q5) * u2 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + 2 * lbo * pxp * (sin(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - sin(q6) * cos(q2) * cos(q3) * cos(q5) - sin(q3) * cos(q2) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - pxp^2 - pyp^2 - pzp^2 - la^2 * u2^2 - la^2 * u3^2 - 2 * la * pyp * sin(q1) * cos(q2) * u2 - 2 * la * pzp * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * u3 - 2 * la * pyp * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * u3 - lbo^2 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)^2 - 2 * la * lbo * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u2 * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - lbo^2 * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)^2 - 2 * lbo * pxp * (sin(q2) * cos(q4) * cos(q5) + sin(q5) * cos(q2) * cos(q3) + sin(q3) * sin(q4) * cos(q2) * cos(q5)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 2 * lbo * pyp * (sin(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) + sin(q6) * cos(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - 2 * lbo * pzp * (cos(q1) * cos(q2) * cos(q4) * cos(q5) + sin(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - sin(q4) * cos(q5) * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1))) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 2 * lbo * pyp * (sin(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - sin(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q4) * cos(q5) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - 2 * lbo * pzp * (sin(q6) * cos(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - cos(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3))

end

function potentialenergy(sol, t)
    @unpack la, lao, lbo, g, fz, ma, mb = sol.prob.p
    q1, q2, q3, q4, q5, q6 = sol(t)

    pz = evaluate(fz, t)

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
    @unpack la, lao, lbo, ixa, ixb, iya, iyb, iza, izb, fx, fy, fz, ma, mb = sol.prob.p
    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = sol(t)

    # pxp, pyp, pzp = sol.prob.p.fp(t)
    pxp = derivative(fx, t, 1)
    pyp = derivative(fy, t, 1)
    pzp = derivative(fz, t, 1)

    amomx = lao * ma * (pzp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - pyp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3))) + mb * (la * pzp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - la * pyp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - lbo * pyp * (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - lbo * pzp * (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)))) + sin(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 + cos(q2) * cos(q3) * (ixa * u1 - la * lbo * mb * ((sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3)) + (sin(q2) * cos(q4) * cos(q5) + sin(q5) * cos(q2) * cos(q3) + sin(q3) * sin(q4) * cos(q2) * cos(q5)) * (izb * u6 + izb * sin(q5) * u1 + izb * cos(q4) * cos(q5) * u3 - izb * sin(q4) * cos(q5) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)) + (sin(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - sin(q6) * cos(q2) * cos(q3) * cos(q5) - sin(q3) * cos(q2) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (iyb * u5 + iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 + iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 - iyb * sin(q6) * cos(q5) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)) + (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) * (ixb * u4 + ixb * cos(q5) * cos(q6) * u1 + ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 + la * lbo * mb * (sin(q5) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - sin(q6) * cos(q5) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2))) - sin(q3) * cos(q2) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2
    amomy = (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 + (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) * (ixa * u1 - la * lbo * mb * ((sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3)) + (sin(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - sin(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q4) * cos(q5) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (izb * u6 + izb * sin(q5) * u1 + izb * cos(q4) * cos(q5) * u3 - izb * sin(q4) * cos(q5) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)) - lao * ma * (pzp * cos(q2) * cos(q3) - pxp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3))) - mb * (la * pzp * cos(q2) * cos(q3) + lbo * pzp * (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - la * pxp * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - lbo * pxp * (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)))) - sin(q1) * cos(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 - (sin(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) + sin(q6) * cos(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (iyb * u5 + iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 + iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 - iyb * sin(q6) * cos(q5) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)) - (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (ixb * u4 + ixb * cos(q5) * cos(q6) * u1 + ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 + la * lbo * mb * (sin(q5) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - sin(q6) * cos(q5) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2)))
    amomz = lao * ma * (pyp * cos(q2) * cos(q3) - pxp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3))) + cos(q1) * cos(q2) * (iza + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u3 + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (iya + ma * lao^2 + la * mb * (la + lbo * cos(q5) * cos(q6))) * u2 + (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) * (ixa * u1 - la * lbo * mb * ((sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3)) + (cos(q1) * cos(q2) * cos(q4) * cos(q5) + sin(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - sin(q4) * cos(q5) * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1))) * (izb * u6 + izb * sin(q5) * u1 + izb * cos(q4) * cos(q5) * u3 - izb * sin(q4) * cos(q5) * u2 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3)) + (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) * (ixb * u4 + ixb * cos(q5) * cos(q6) * u1 + ixb * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * u2 + ixb * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) * u3 + la * lbo * mb * (sin(q5) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - sin(q6) * cos(q5) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2))) - mb * (la * pxp * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - la * pyp * cos(q2) * cos(q3) - lbo * pyp * (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - lbo * pxp * (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)))) - (sin(q6) * cos(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - cos(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (iyb * u5 + iyb * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 + iyb * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2 - iyb * sin(q6) * cos(q5) * u1 - lbo * mb * (lbo + la * cos(q5) * cos(q6)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2))
   
    return amomx, amomy, amomz
end

angmom(sol) = [tup[k] for tup ∈ [angmom(sol, t) for t ∈ sol.t], k ∈ 1:3]
