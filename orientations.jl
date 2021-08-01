# Calculates orientation angles, derivatives and generalised speeds
# Body123
using Rotations

function getorientation(time, base, mid, tip)
    # Set up vectors
    v1 = normalize(mid - base, 2); v2 = normalize(tip - mid, 2)

    # Get rotation matricies
    xaxis = repeat([1 0 0], size(base, 1))
    r1 = [(rotation_between(xaxis[i,:], v1[i,:]) |> RotXYZ) for i in 1:size(v1, 1)]
    r2 = [(rotation_between(v1[i,:], v2[i,:]) |> RotXYZ) for i in 1:size(v2, 1)]

    # Extract Euler angles
    q1 = getfield.(r1, :theta1); q2 = getfield.(r1, :theta2); q3 = getfield.(r1, :theta3)
    q4 = getfield.(r2, :theta1); q5 = getfield.(r2, :theta2); q6 = getfield.(r2, :theta3)

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
    
    return [q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6]

end

# Get XYZ intrinsic Euler angles
function getorientation(v::Vector{Float64}, ref=[1,0,0])
    # Create vector
    v_norm = normalize(v)

    # Create rotation matrix
    r = rotation_between(ref, v_norm) |> RotXYZ

    # Extract angles
    θ₁, θ₂, θ₃ = r.theta1, r.theta2, r.theta3 

    return [θ₁, θ₂, θ₃]
end

function getorientation(m::Matrix{Float64}, ref::Vector{Vector{Float64}})
    θ = similar(m)
    @inbounds for i ∈ 1:size(m, 1)
        θ[i,:] = getorientation(m[i,:], ref[i])
    end

    return θ
end