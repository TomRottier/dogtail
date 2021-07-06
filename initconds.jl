# Calculates initial conditions; numerically solves inverse kinematics
using NLsolve

function getinitcond(initconds, p)

    px, py, pz, pxp, pyp, pzp, p2x, p2y, p2z, vp2x, vp2y, vp2z, p3x, p3y, p3z, vp3x, vp3y, vp3z = initconds
    la = p.la; lb = p.lb
    
    function f!(F, x)
        q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = x
        
        # F = Vector{Float64}(undef, 12)
        F[1] = px + la * cos(q2) * cos(q3) - p2x
        F[2] = py + la * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - p2y
        F[3] = pz + la * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - p2z
        F[4] = px + la * cos(q2) * cos(q3) + lb * (cos(q2) * cos(q3) * cos(q5) * cos(q6) + sin(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - sin(q3) * cos(q2) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - p3x
        F[5] = py + la * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - lb * (sin(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) - cos(q5) * cos(q6) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6)) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) - p3y
        F[6] = pz + la * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + lb * (cos(q1) * cos(q2) * (sin(q4) * sin(q6) - sin(q5) * cos(q4) * cos(q6)) + cos(q5) * cos(q6) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) + (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (sin(q6) * cos(q4) + sin(q4) * sin(q5) * cos(q6))) - p3z
        F[7] = pxp - la * sin(q2) * u2 - la * sin(q3) * cos(q2) * u3 - vp2x
        F[8] = pyp + la * sin(q1) * cos(q2) * u2 + la * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * u3 - vp2y
        F[9] = pzp + la * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * u3 - la * cos(q1) * cos(q2) * u2 - vp2z
        F[10] = pxp + lb * (sin(q2) * cos(q4) * cos(q5) + sin(q5) * cos(q2) * cos(q3) + sin(q3) * sin(q4) * cos(q2) * cos(q5)) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - la * sin(q2) * u2 - la * sin(q3) * cos(q2) * u3 - lb * (sin(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - sin(q6) * cos(q2) * cos(q3) * cos(q5) - sin(q3) * cos(q2) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - vp3x
        F[11] = pyp + la * sin(q1) * cos(q2) * u2 + la * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * u3 + lb * (sin(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) + sin(q6) * cos(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) + lb * (sin(q5) * (sin(q3) * cos(q1) + sin(q1) * sin(q2) * cos(q3)) - sin(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q4) * cos(q5) * (cos(q1) * cos(q3) - sin(q1) * sin(q2) * sin(q3))) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) - vp3y
        F[12] = pzp + la * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * u3 + lb * (cos(q1) * cos(q2) * cos(q4) * cos(q5) + sin(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - sin(q4) * cos(q5) * (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1))) * (sin(q6) * cos(q5) * u1 - u5 - (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) * u3 - (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6)) * u2) + lb * (sin(q6) * cos(q5) * (sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)) - cos(q1) * cos(q2) * (sin(q4) * cos(q6) + sin(q5) * sin(q6) * cos(q4)) - (sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1)) * (cos(q4) * cos(q6) - sin(q4) * sin(q5) * sin(q6))) * (sin(q4) * cos(q5) * u2 - u6 - sin(q5) * u1 - cos(q4) * cos(q5) * u3) - la * cos(q1) * cos(q2) * u2 - vp3z
      end

    # x₀ = repeat([0.], 12)
    # options = Options(func=f!, N=12, Ns=50, Nt=200, lb=repeat([-7.], 12), ub=repeat([7.], 12), print_status=false, c=repeat([2.0], 12))
    # result = Result(fopt=options.f(x₀), xopt=x₀)
    # current = State(f=options.f(x₀), x=x₀, v=options.ub .- options.lb, T=100.0)

    # sa!(current, result, options)



    sol = nlsolve(f!, repeat([0.], 12), method=:trust_region, store_trace=true, extended_trace=true, autodiff=:forward)

    q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6 = sol.zero
    
    return [q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6]
end


# space123
# F[1]  = px + la * cos(q2) * cos(q3) - p2x
# F[2]  = py + la * sin(q3) * cos(q2) - p2y
# F[3]  = pz - la * sin(q2) - p2z
# F[4]  = px + la * cos(q2) * cos(q3) + lb * (cos(q2) * cos(q3) * cos(q5) * cos(q6) - sin(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - sin(q6) * cos(q5) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3))) - p3x
# F[5]  = py + la * sin(q3) * cos(q2) + lb * (sin(q3) * cos(q2) * cos(q5) * cos(q6) + sin(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) + sin(q6) * cos(q5) * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3))) - p3y
# F[6]  = pz - la * sin(q2) - lb * (sin(q2) * cos(q5) * cos(q6) + sin(q5) * cos(q1) * cos(q2) - sin(q1) * sin(q6) * cos(q2) * cos(q5)) - p3z
# F[7]  = pxp - la * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) * u2 - la * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * u3 - vp2x
# F[8]  = pyp + la * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * u3 + la * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) * u2 - vp2y
# F[9]  = pzp + la * sin(q1) * cos(q2) * u3 - la * cos(q1) * cos(q2) * u2 - vp2z
# F[10] = pxp + lb * (sin(q4) * cos(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) - cos(q2) * cos(q3) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) - (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3))) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - la * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) * u2 - la * (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * u3 - lb * (cos(q2) * cos(q3) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) + cos(q4) * cos(q5) * (sin(q1) * sin(q3) + sin(q2) * cos(q1) * cos(q3)) + (sin(q3) * cos(q1) - sin(q1) * sin(q2) * cos(q3)) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) - vp3x
# F[11] = pyp + la * (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * u3 + la * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) * u2 - lb * (sin(q3) * cos(q2) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) + sin(q4) * cos(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) - (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6))) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - lb * (sin(q3) * cos(q2) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) - cos(q4) * cos(q5) * (sin(q1) * cos(q3) - sin(q2) * sin(q3) * cos(q1)) - (cos(q1) * cos(q3) + sin(q1) * sin(q2) * sin(q3)) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) - vp3y
# F[12] = pzp + la * sin(q1) * cos(q2) * u3 + lb * (sin(q4) * cos(q1) * cos(q2) * cos(q5) + sin(q2) * (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) + sin(q1) * cos(q2) * (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6))) * (u6 + cos(q4) * cos(q5) * u3 + (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) * u1 - (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4)) * u2) - la * cos(q1) * cos(q2) * u2 - lb * (cos(q1) * cos(q2) * cos(q4) * cos(q5) - sin(q2) * (sin(q4) * sin(q6) + sin(q5) * cos(q4) * cos(q6)) - sin(q1) * cos(q2) * (sin(q4) * cos(q6) - sin(q5) * sin(q6) * cos(q4))) * (u5 + sin(q4) * cos(q5) * u3 + (cos(q4) * cos(q6) + sin(q4) * sin(q5) * sin(q6)) * u2 - (sin(q6) * cos(q4) - sin(q4) * sin(q5) * cos(q6)) * u1) - vp3z


# body123
# F[1] = px + la*cos(q2)*cos(q3) - p2x
# F[2] = py + la*(sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3)) - p2y
# F[3] = pz + la*(sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3)) - p2z
# F[4] = px + la*cos(q2)*cos(q3) + lb*(cos(q2)*cos(q3)*cos(q5)*cos(q6)+sin(q2)*(sin(q4)*sin(q6)-sin(q5)*cos(q4)*cos(q6))-sin(q3)*cos(q2)*(sin(q6)*cos(q4)+sin(q4)*sin(q5)*cos(q6))) - p3x
# F[5] = py + la*(sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3)) - lb*(sin(q1)*cos(q2)*(sin(q4)*sin(q6)-sin(q5)*cos(q4)*cos(q6))-cos(q5)*cos(q6)*(sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))-(sin(q6)*cos(q4)+sin(q4)*sin(q5)*cos(q6))*(cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))) - p3y
# F[6] = pz + la*(sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3)) + lb*(cos(q1)*cos(q2)*(sin(q4)*sin(q6)-sin(q5)*cos(q4)*cos(q6))+cos(q5)*cos(q6)*(sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))+(sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*(sin(q6)*cos(q4)+sin(q4)*sin(q5)*cos(q6))) - p3z
# F[7] = pxp - la*sin(q2)*u2 - la*sin(q3)*cos(q2)*u3 - vp2x
# F[8] = pyp + la*sin(q1)*cos(q2)*u2 + la*(cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*u3 - vp2y
# F[9] = pzp + la*(sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*u3 - la*cos(q1)*cos(q2)*u2 - vp2z
# F[10] = pxp + lb*(sin(q2)*cos(q4)*cos(q5)+sin(q5)*cos(q2)*cos(q3)+sin(q3)*sin(q4)*cos(q2)*cos(q5))*(sin(q6)*cos(q5)*u1-u5-(sin(q4)*cos(q6)+sin(q5)*sin(q6)*cos(q4))*u3-(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))*u2) - la*sin(q2)*u2 - la*sin(q3)*cos(q2)*u3 - lb*(sin(q2)*(sin(q4)*cos(q6)+sin(q5)*sin(q6)*cos(q4))-sin(q6)*cos(q2)*cos(q3)*cos(q5)-sin(q3)*cos(q2)*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6)))*(sin(q4)*cos(q5)*u2-u6-sin(q5)*u1-cos(q4)*cos(q5)*u3) - vp3x
# F[11] = pyp + la*sin(q1)*cos(q2)*u2 + la*(cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*u3 + lb*(sin(q1)*cos(q2)*(sin(q4)*cos(q6)+sin(q5)*sin(q6)*cos(q4))+sin(q6)*cos(q5)*(sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))-(cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6)))*(sin(q4)*cos(q5)*u2-u6-sin(q5)*u1-cos(q4)*cos(q5)*u3) + lb*(sin(q5)*(sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))-sin(q1)*cos(q2)*cos(q4)*cos(q5)-sin(q4)*cos(q5)*(cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3)))*(sin(q6)*cos(q5)*u1-u5-(sin(q4)*cos(q6)+sin(q5)*sin(q6)*cos(q4))*u3-(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))*u2) - vp3y
# F[12] = pzp + la*(sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*u3 + lb*(cos(q1)*cos(q2)*cos(q4)*cos(q5)+sin(q5)*(sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))-sin(q4)*cos(q5)*(sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1)))*(sin(q6)*cos(q5)*u1-u5-(sin(q4)*cos(q6)+sin(q5)*sin(q6)*cos(q4))*u3-(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))*u2) + lb*(sin(q6)*cos(q5)*(sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))-cos(q1)*cos(q2)*(sin(q4)*cos(q6)+sin(q5)*sin(q6)*cos(q4))-(sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6)))*(sin(q4)*cos(q5)*u2-u6-sin(q5)*u1-cos(q4)*cos(q5)*u3) - la*cos(q1)*cos(q2)*u2 - vp3z
