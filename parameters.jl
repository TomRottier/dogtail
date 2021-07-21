# Parameters struct
@with_kw mutable struct Params{T} @deftype T
    la
    lb
    lao
    lbo
    ma = 0.3
    mb = 0.3
    r = 0.01
    ixa = 0.5 * ma * r^2
    ixb = 0.5 * mb * r^2
    iya = 1 / 12 * ma * (3 * r^2 + la^2)
    iyb = 1 / 12 * mb * (3 * r^2 + lb^2)
    iza = 1 / 12 * ma * (3 * r^2 + la^2)
    izb = 1 / 12 * mb * (3 * r^2 + lb^2)
    g = -9.81
    fx::Spline1D
    fy::Spline1D
    fz::Spline1D
    ka = 0.01
    kb = 0.01
    ba = 0.001
    bb = 0.001
    eqx = 0.0
    eqy = 0.0
    eqz = 0.0
end

# Set parameters for model
function initialise_parameters(fx, fy, fz, la, lb)
    # Calculate inertial parameters of tail segments assuming solid cylinders of uniform density
    lao = 0.5 * la; lbo = 0.5 * lb

    return Params(la=la, lb=lb, lao=lao, lbo=lbo, fx=fx, fy=fy, fz=fz)
end

function update_parameters!(pin, p::Params)
    @inbounds ma, mb, ka, kb, ba, bb, eqx, eqy, eqz = pin
    @unpack la, lb, r = p
    # Update inertial parameters of tail segments from masses
    ixa = 0.5 * ma * r^2; ixb = 0.5 * mb * r^2
    iya = iza = 1 / 12 * ma * (3 * r^2 + la^2)
    iyb = izb = 1 / 12 * mb * (3 * r^2 + lb^2)

    return  @pack! p = ixa, ixb, iya, iyb, iza, izb, ma, mb, ka, kb, ba, bb, eqx, eqy, eqz
end


