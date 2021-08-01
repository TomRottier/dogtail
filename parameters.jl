# File generated automatically from dogtail.f
# Parameters struct
@with_kw mutable struct Params{T} @deftype T
    z::Vector{Float64} = Vector{Float64}(undef, 268)
    la
    lb
    lao = 0.5 * la
    lbo = 0.5 * lb
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
    ox::Union{Float64,Spline1D}
	oy::Union{Float64,Spline1D}
	oz::Union{Float64,Spline1D}
	px::Union{Float64,Spline1D}
	py::Union{Float64,Spline1D}
	pz::Union{Float64,Spline1D}
    ka = 0.01
    kb = 0.01
    ba = 0.001
    bb = 0.001
    eqx = 0.0
    eqy = 0.0
    eqz = 0.0
end

# Set parameters for model
function initialise_parameters(px, py, pz, ox, oy, oz, la, lb)
    # Calculate inertial parameters of tail segments assuming solid cylinders of uniform density
    lao = 0.5 * la; lbo = 0.5 * lb
    p = Params(la=la, lb=lb, lao=lao, lbo=lbo, px=px, py=py, pz=pz, ox=ox, oy=oy, oz=oz)
    @unpack g, ixb, iya, iyb, iza, izb, la, lao, lb, lbo, ma, mb, z = p

    # Evaluate constants
    u4 = 0
	u4p = 0
	z[135] = g * ma
	z[136] = g * mb
	z[156] = lao * z[135]
	z[159] = lbo * z[136]
	z[180] = ixb * u4
	z[219] = lbo * mb
	z[223] = iya + ma * lao^2
	z[224] = la^2
	z[230] = lao * ma
	z[232] = iza + ma * lao^2
	z[239] = iyb + mb * lbo^2
	z[244] = izb + mb * lbo^2
	z[255] = ma + mb
	z[260] = la * mb

    return p
end

function update_parameters!(p::Params, pin)
    @inbounds ma, mb, ka, kb, ba, bb, eqx, eqy, eqz = pin
    @unpack la, lb, r, g, lao, lbo, la,  z = p

    # Update inertial parameters of tail segments from masses
    ixa = 0.5 * ma * r^2; ixb = 0.5 * mb * r^2
    iya = iza = 1 / 12 * ma * (3 * r^2 + la^2)
    iyb = izb = 1 / 12 * mb * (3 * r^2 + lb^2)

    # Update constants
    u4 = 0
	u4p = 0
	z[135] = g * ma
	z[136] = g * mb
	z[156] = lao * z[135]
	z[159] = lbo * z[136]
	z[180] = ixb * u4
	z[219] = lbo * mb
	z[223] = iya + ma * lao^2
	z[224] = la^2
	z[230] = lao * ma
	z[232] = iza + ma * lao^2
	z[239] = iyb + mb * lbo^2
	z[244] = izb + mb * lbo^2
	z[255] = ma + mb
	z[260] = la * mb

    return  @pack! p = ixa, ixb, iya, iyb, iza, izb, ma, mb, ka, kb, ba, bb, eqx, eqy, eqz
end
