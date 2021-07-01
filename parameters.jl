# Set parameters for model
function getparameters(pin)
    # Functions and lengths fixed for given data; stiffness, damping and masses change
    fx, fy, fz, la, lb, ka, kb, ba, bb, ma, mb = pin

    # Calculate other variables
    lao = 0.5 * la; lbo = 0.5 * lb
    ixa = ixb = 0.01;
    iya = iza = 1 / 12 * ma * lao;
    iyb = izb = 1 / 12 * mb * lbo
    g = -9.81

    return (la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, fx = fx, fy = fy, fz = fz, ka, kb, ba, bb, ma, mb)
end
