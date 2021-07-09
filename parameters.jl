# Set parameters for model
function initialise_parameters(pinitial)
    # Functions and lengths fixed for given data; stiffness, damping and masses change
    fx, fy, fz, la, lb, ka, kb, ba, bb, ma, mb = pinitial

    # Calculate inertial parameters of tail segments assuming solid cylinders of uniform density
    lao = 0.5 * la; lbo = 0.5 * lb
    r = 0.01    # 2 cm diameter tail
    ixa = 0.5 * ma * r^2; ixb = 0.5 * mb * r^2;
    iya = iza = 1 / 12 * ma * (3 * r^2 + la^2);
    iyb = izb = 1 / 12 * mb * (3 * r^2 + lb^2)

    # Fixed parameters
    g = -9.81

    return (la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx, fy, fz, ka, kb, ba, bb, ma, mb)
end

function update_parameters(p)
    # Update inertial parameters of tail segments from masses
    # ma,mb = p
    r = 0.01    # 2 cm diameter tail
    ixa = 0.5 * ma * r^2; ixb = 0.5 * mb * r^2;
    iya = iza = 1 / 12 * ma * (3 * r^2 + la^2);
    iyb = izb = 1 / 12 * mb * (3 * r^2 + lb^2)

    return (ixa, ixb, iya, iyb, iza, izb)
end


