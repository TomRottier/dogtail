# Equations of motion for 3d pendulum
function eom(du, u, p, t)
    q4, q5, q6, u4, u5, u6 = u
    
    q4p = u4 + tan(q5) * (sin(q4) * u5 + cos(q4) * u6)
    q5p = cos(q4) * u5 - sin(q4) * u6
    q6p = (sin(q4) * u5 + cos(q4) * u6) / cos(q5)

    du[1] = q4p
    du[2] = q5p
    du[3] = q6p
    du[4:6] = acc(u, p)

end


function acc(u, p)
    ma, ix, iy, iz, lao, la, g = p

    q4, q5, q6, u4, u5, u6 = u

    coef = Matrix{Float64}(undef, 3, 3)
    rhs = Vector{Float64}(undef, 3)

    coef[1,1] = -ix
    coef[1,2] = 0
    coef[1,3] = 0
    coef[2,1] = 0
    coef[2,2] = -iy - m * lao^2
    coef[2,3] = 0
    coef[3,1] = 0
    coef[3,2] = 0
    coef[3,3] = -iz - m * lao^2

    rhs[1] = -(iy - iz) * u5 * u6
    rhs[2] = g * lao * m * cos(q4) * cos(q5) + (ix - iz - m * lao^2) * u4 * u6
    rhs[3] = -g * lao * m * sin(q4) * cos(q5) - (ix - iy - m * lao^2) * u4 * u5
    
    return coef \ rhs

end


# kindiffs generated
# q4p = u4 + tan(q5)*(sin(q4)*u5+cos(q4)*u6)
# q5p = cos(q4)*u5 - sin(q4)*u6
# q6p = (sin(q4)*u5+cos(q4)*u6)/cos(q5)

# coef[1,1] = -ix
# coef[1,2] = 0
# coef[1,3] = 0
# coef[2,1] = 0
# coef[2,2] = -iy - m*lao^2
# coef[2,3] = 0
# coef[3,1] = 0
# coef[3,2] = 0
# coef[3,3] = -iz - m*lao^2
# rhs[1] = -(iy-iz)*u5*u6
# rhs[2] = g*lao*m*cos(q4)*cos(q5) + (ix-iz-m*lao^2)*u4*u6
# rhs[3] = -g*lao*m*sin(q4)*cos(q5) - (ix-iy-m*lao^2)*u4*u5


# angvel expressed in n axes
# q4p = sin(q6)*cos(q5)*u5 + cos(q5)*cos(q6)*u4 + tan(q5)*(sin(q4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)+cos(q4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5)) - sin(q5)*u6
# q5p = cos(q4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4) - sin(q4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5)
# q6p = (sin(q4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)+cos(q4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5))/cos(q5)

# coef[1,1] = -ix*cos(q5)^2*cos(q6)^2 - iz*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))^2 - iy*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))^2 - m*lao^2*((sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))^2+(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))^2)
# coef[1,2] = iy*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6)) + iz*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4)) + m*lao^2*((sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))) - ix*sin(q6)*cos(q6)*cos(q5)^2
# coef[1,3] = cos(q5)*(ix*sin(q5)*cos(q6)+iy*sin(q4)*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))-iz*cos(q4)*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))-m*lao^2*(cos(q4)*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))-sin(q4)*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))))
# coef[2,1] = iy*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6)) + iz*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4)) + m*lao^2*((sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))) - ix*sin(q6)*cos(q6)*cos(q5)^2
# coef[2,2] = -ix*sin(q6)^2*cos(q5)^2 - iy*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))^2 - iz*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))^2 - m*lao^2*((cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))^2+(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))^2)
# coef[2,3] = cos(q5)*(ix*sin(q5)*sin(q6)+iz*cos(q4)*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))-iy*sin(q4)*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))-m*lao^2*(sin(q4)*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))-cos(q4)*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))))
# coef[3,1] = cos(q5)*(ix*sin(q5)*cos(q6)+iy*sin(q4)*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))-iz*cos(q4)*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))-m*lao^2*(cos(q4)*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))-sin(q4)*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))))
# coef[3,2] = cos(q5)*(ix*sin(q5)*sin(q6)+iz*cos(q4)*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))-iy*sin(q4)*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))-m*lao^2*(sin(q4)*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))-cos(q4)*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))))
# coef[3,3] = -ix*sin(q5)^2 - m*lao^2*cos(q5)^2 - iy*sin(q4)^2*cos(q5)^2 - iz*cos(q4)^2*cos(q5)^2
# rhs[1] = (ix-iy)*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*(sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4) + (ix-iz)*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*(sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5) - g*lao*m*cos(q5)*(sin(q4)*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))+cos(q4)*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))) - (iy-iz)*cos(q5)*cos(q6)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5) - ix*cos(q5)*cos(q6)*(cos(q5)*u6*u5+sin(q5)*sin(q6)*u5*u5+sin(q5)*cos(q6)*u4*u5+sin(q6)*cos(q5)*u4*u6-cos(q5)*cos(q6)*u5*u6) - iz*(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*(sin(q4)*cos(q5)*u6*u4+sin(q5)*cos(q4)*u6*u5+u4*(sin(q4)*sin(q5)*cos(q6)*u4+sin(q5)*sin(q6)*cos(q4)*u6-sin(q4)*cos(q6)*u6-sin(q6)*cos(q4)*u4-cos(q4)*cos(q5)*cos(q6)*u5)+u5*(cos(q4)*cos(q6)*u4+sin(q4)*sin(q5)*sin(q6)*u4-sin(q4)*sin(q6)*u6-sin(q5)*cos(q4)*cos(q6)*u6-sin(q6)*cos(q4)*cos(q5)*u5)) - iy*(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*(cos(q4)*cos(q5)*u6*u4-sin(q4)*sin(q5)*u6*u5-u4*(cos(q4)*cos(q6)*u6+sin(q4)*sin(q5)*sin(q6)*u6-sin(q4)*sin(q6)*u4-sin(q4)*cos(q5)*cos(q6)*u5-sin(q5)*cos(q4)*cos(q6)*u4)-u5*(sin(q4)*cos(q6)*u4+sin(q6)*cos(q4)*u6-sin(q4)*sin(q5)*cos(q6)*u6-sin(q4)*sin(q6)*cos(q5)*u5-sin(q5)*sin(q6)*cos(q4)*u4)) - m*lao^2*((sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*((sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)+sin(q4)*cos(q5)*u6*u4+sin(q5)*cos(q4)*u6*u5+u4*(sin(q4)*sin(q5)*cos(q6)*u4+sin(q5)*sin(q6)*cos(q4)*u6-sin(q4)*cos(q6)*u6-sin(q6)*cos(q4)*u4-cos(q4)*cos(q5)*cos(q6)*u5)+u5*(cos(q4)*cos(q6)*u4+sin(q4)*sin(q5)*sin(q6)*u4-sin(q4)*sin(q6)*u6-sin(q5)*cos(q4)*cos(q6)*u6-sin(q6)*cos(q4)*cos(q5)*u5))+(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*((sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5)+cos(q4)*cos(q5)*u6*u4-sin(q4)*sin(q5)*u6*u5-u4*(cos(q4)*cos(q6)*u6+sin(q4)*sin(q5)*sin(q6)*u6-sin(q4)*sin(q6)*u4-sin(q4)*cos(q5)*cos(q6)*u5-sin(q5)*cos(q4)*cos(q6)*u4)-u5*(sin(q4)*cos(q6)*u4+sin(q6)*cos(q4)*u6-sin(q4)*sin(q5)*cos(q6)*u6-sin(q4)*sin(q6)*cos(q5)*u5-sin(q5)*sin(q6)*cos(q4)*u4)))
# rhs[2] = g*lao*m*cos(q5)*(cos(q4)*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))+sin(q4)*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))) + iz*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*(sin(q4)*cos(q5)*u6*u4+sin(q5)*cos(q4)*u6*u5+u4*(sin(q4)*sin(q5)*cos(q6)*u4+sin(q5)*sin(q6)*cos(q4)*u6-sin(q4)*cos(q6)*u6-sin(q6)*cos(q4)*u4-cos(q4)*cos(q5)*cos(q6)*u5)+u5*(cos(q4)*cos(q6)*u4+sin(q4)*sin(q5)*sin(q6)*u4-sin(q4)*sin(q6)*u6-sin(q5)*cos(q4)*cos(q6)*u6-sin(q6)*cos(q4)*cos(q5)*u5)) + iy*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*(cos(q4)*cos(q5)*u6*u4-sin(q4)*sin(q5)*u6*u5-u4*(cos(q4)*cos(q6)*u6+sin(q4)*sin(q5)*sin(q6)*u6-sin(q4)*sin(q6)*u4-sin(q4)*cos(q5)*cos(q6)*u5-sin(q5)*cos(q4)*cos(q6)*u4)-u5*(sin(q4)*cos(q6)*u4+sin(q6)*cos(q4)*u6-sin(q4)*sin(q5)*cos(q6)*u6-sin(q4)*sin(q6)*cos(q5)*u5-sin(q5)*sin(q6)*cos(q4)*u4)) + m*lao^2*((sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*((sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)+sin(q4)*cos(q5)*u6*u4+sin(q5)*cos(q4)*u6*u5+u4*(sin(q4)*sin(q5)*cos(q6)*u4+sin(q5)*sin(q6)*cos(q4)*u6-sin(q4)*cos(q6)*u6-sin(q6)*cos(q4)*u4-cos(q4)*cos(q5)*cos(q6)*u5)+u5*(cos(q4)*cos(q6)*u4+sin(q4)*sin(q5)*sin(q6)*u4-sin(q4)*sin(q6)*u6-sin(q5)*cos(q4)*cos(q6)*u6-sin(q6)*cos(q4)*cos(q5)*u5))+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*((sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5)+cos(q4)*cos(q5)*u6*u4-sin(q4)*sin(q5)*u6*u5-u4*(cos(q4)*cos(q6)*u6+sin(q4)*sin(q5)*sin(q6)*u6-sin(q4)*sin(q6)*u4-sin(q4)*cos(q5)*cos(q6)*u5-sin(q5)*cos(q4)*cos(q6)*u4)-u5*(sin(q4)*cos(q6)*u4+sin(q6)*cos(q4)*u6-sin(q4)*sin(q5)*cos(q6)*u6-sin(q4)*sin(q6)*cos(q5)*u5-sin(q5)*sin(q6)*cos(q4)*u4))) - (ix-iz)*(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*(sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5) - (ix-iy)*(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*(sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4) - (iy-iz)*sin(q6)*cos(q5)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5) - ix*sin(q6)*cos(q5)*(cos(q5)*u6*u5+sin(q5)*sin(q6)*u5*u5+sin(q5)*cos(q6)*u4*u5+sin(q6)*cos(q5)*u4*u6-cos(q5)*cos(q6)*u5*u6)
# rhs[3] = (ix-iy)*cos(q4)*cos(q5)*(sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4) + (iy-iz)*sin(q5)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5) + ix*sin(q5)*(cos(q5)*u6*u5+sin(q5)*sin(q6)*u5*u5+sin(q5)*cos(q6)*u4*u5+sin(q6)*cos(q5)*u4*u6-cos(q5)*cos(q6)*u5*u6) + iy*sin(q4)*cos(q5)*(cos(q4)*cos(q5)*u6*u4-sin(q4)*sin(q5)*u6*u5-u4*(cos(q4)*cos(q6)*u6+sin(q4)*sin(q5)*sin(q6)*u6-sin(q4)*sin(q6)*u4-sin(q4)*cos(q5)*cos(q6)*u5-sin(q5)*cos(q4)*cos(q6)*u4)-u5*(sin(q4)*cos(q6)*u4+sin(q6)*cos(q4)*u6-sin(q4)*sin(q5)*cos(q6)*u6-sin(q4)*sin(q6)*cos(q5)*u5-sin(q5)*sin(q6)*cos(q4)*u4)) - (ix-iz)*sin(q4)*cos(q5)*(sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5) - iz*cos(q4)*cos(q5)*(sin(q4)*cos(q5)*u6*u4+sin(q5)*cos(q4)*u6*u5+u4*(sin(q4)*sin(q5)*cos(q6)*u4+sin(q5)*sin(q6)*cos(q4)*u6-sin(q4)*cos(q6)*u6-sin(q6)*cos(q4)*u4-cos(q4)*cos(q5)*cos(q6)*u5)+u5*(cos(q4)*cos(q6)*u4+sin(q4)*sin(q5)*sin(q6)*u4-sin(q4)*sin(q6)*u6-sin(q5)*cos(q4)*cos(q6)*u6-sin(q6)*cos(q4)*cos(q5)*u5)) - m*lao^2*cos(q5)*(cos(q4)*((sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(sin(q4)*cos(q5)*u6+(cos(q4)*cos(q6)+sin(q4)*sin(q5)*sin(q6))*u5-(sin(q6)*cos(q4)-sin(q4)*sin(q5)*cos(q6))*u4)+sin(q4)*cos(q5)*u6*u4+sin(q5)*cos(q4)*u6*u5+u4*(sin(q4)*sin(q5)*cos(q6)*u4+sin(q5)*sin(q6)*cos(q4)*u6-sin(q4)*cos(q6)*u6-sin(q6)*cos(q4)*u4-cos(q4)*cos(q5)*cos(q6)*u5)+u5*(cos(q4)*cos(q6)*u4+sin(q4)*sin(q5)*sin(q6)*u4-sin(q4)*sin(q6)*u6-sin(q5)*cos(q4)*cos(q6)*u6-sin(q6)*cos(q4)*cos(q5)*u5))-sin(q4)*((sin(q5)*u6-sin(q6)*cos(q5)*u5-cos(q5)*cos(q6)*u4)*(cos(q4)*cos(q5)*u6+(sin(q4)*sin(q6)+sin(q5)*cos(q4)*cos(q6))*u4-(sin(q4)*cos(q6)-sin(q5)*sin(q6)*cos(q4))*u5)+cos(q4)*cos(q5)*u6*u4-sin(q4)*sin(q5)*u6*u5-u4*(cos(q4)*cos(q6)*u6+sin(q4)*sin(q5)*sin(q6)*u6-sin(q4)*sin(q6)*u4-sin(q4)*cos(q5)*cos(q6)*u5-sin(q5)*cos(q4)*cos(q6)*u4)-u5*(sin(q4)*cos(q6)*u4+sin(q6)*cos(q4)*u6-sin(q4)*sin(q5)*cos(q6)*u6-sin(q4)*sin(q6)*cos(q5)*u5-sin(q5)*sin(q6)*cos(q4)*u4)))
