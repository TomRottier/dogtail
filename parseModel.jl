## Read whole file as string
fname = "dogtail.f"
lines = read("model/$fname", String)

# Convert string from Fortran to Julia
function convertToJulia(x::AbstractString)
    # Remove continuation characters
    out = replace(x, r"\n\s+&" => "")

    # Remove whitespace
    out = replace(out, r"^\s+" => "") |> x -> replace(x, r"\r" => "") |> x -> replace(x, r"\n" => "")

    # Change index notation
    out = replace(out, r"\s*(\w+)\((([0-9]){1,3})\) = " => s"\1[\2] = ")
    # Change matrix index notation
    out = replace(out, r"\s*(\w+)\((([0-9]){1,3},([0-9]){1,3})\) = " => s"\1[\2] = ")

    # Change exponential notation
    out = replace(out, "**" => "^")

    # Remove d0's
    out = replace(out, "D0" => "")

    # Make lowercase
    out = lowercase.(out)

    return out
end

# Replace all Z indexes with square brackets
Zindex(x::AbstractString) = replace(x, r"(z|Z)\((([0-9]){1,3})\)" => s"\1[\2]")

## Extract info
# Variables
r = r"COMMON/VARIBLES/ (.*(\n\s+&.*){0,100})"
variables = match(r, lines).captures[1] |>
    x -> split(x, ',') |>
    x -> convertToJulia.(x)

# Parameters
r = r"COMMON/CONSTNTS/ (.*(\n\s+&.*){0,100})"
parameters = match(r, lines).captures[1] |>
    x -> split(x, ',') |>
    x -> convertToJulia.(x)

# Specified
r = r"COMMON/SPECFIED/ (.*(\n\s+&.*){0,100})"
specified = match(r, lines).captures[1] |>
    x -> split(x, ',') |>
    x -> convertToJulia.(x)

# Get the 0th derivative of specified variables as they are used for function names
# All variables not ending in p
specified_functions = specified[.!occursin.(r"p$", specified)]

# Z array
r = r"COMMON/MISCLLNS/ .*(Z\((([0-9]){0,3})\))"
z_array = match(r, lines).captures[2] |> convertToJulia

# Evaluate constants - Z assignments that are constant throughout
r = r"C\*\*\s*Evaluate constants\r\n((.*\n){0,100})C\*\*\s+Initialize"
constants = match(r, lines).captures[1] |>
    x -> Zindex.(x) |>
    x -> split(x, r"\r\n") |>
    x -> convertToJulia.(x) |>
    x -> x[.!isempty.(x)]

# EQNS
rEQNS = r"\s+SUBROUTINE \s* EQNS1"; rIO = r"\s+SUBROUTINE \s* IO"
eqns = lines[match(rEQNS, lines).offset:match(rIO, lines).offset]

# Z and Qp assignments - keeps in correct order
r = r"\s*(Z\(([0-9]{1,3})\)|Q[0-9]p) = .*(\n\s+&.*){0,100}"
z_qp = String[]
for m ∈ eachmatch(r, eqns)
    out = m.match |> convertToJulia |> Zindex
    push!(z_qp, out)
end

# COEF
r = r"\s+COEF.*"
coef = String[]
for m ∈ eachmatch(r, eqns)
    out = m.match |> convertToJulia |> Zindex
    push!(coef, out)
end

# RHS
r = r"\s+RHS.*"
rhs = String[]
for m ∈ eachmatch(r, eqns)
    out = m.match |> convertToJulia |> Zindex
    push!(rhs, out)
end

# IO
io = lines[match(rIO, lines).offset:end]

# Z assignments
r = r"\s*Z\(([0-9]{1,3})\) = .*(\n\s+&.*){0,100}"
z_io = String[]
for m ∈ eachmatch(r, io)
    out = m.match |> convertToJulia |> Zindex
    push!(z_io, out)
end

# All other assignments
r = r"\s+[^Z]\w+ = .*(\n\s+&.*){0,100}\n"
io_all = String[]
for m ∈ eachmatch(r, io)
    out = m.match |> convertToJulia |> Zindex
    push!(io_all, out)
end

# Write IO assignments in terms of generalised coordinates, speeds and parameters
assignments = [constants; z_qp; z_io]
io_raw = String[]
replacements = String[]
# Replaces Z index with what its rhs
function replace_rhs(z)
    # Find where z is defined
    # z_idx = match(r"z\[\d+\]", z).match * " = "
    zn = z
    r = r"z\[\d+\]"
    while any(occursin.(r, zn))
        ms = eachmatch(r"z\[\d+\]", zn)
        for m ∈ ms
            z_idx = m.match * " = "
            idx = occursin.(z_idx, assignments)
            _z = replace(assignments[idx][1], r"z\[([0-9]){1,3}\] = " => "")
            zn = replace(zn, m.match => "($_z)")
        end        
    end    
    return zn
end


# Replace all sin and cos calls with variables, stores replacements
function replaceSinCos(eqn)
    # Match sin/cos calls
    m1 = eachmatch(r"sin\((\w+)\)", eqn) |> collect
    m2 = eachmatch(r"cos\((\w+)\)", eqn) |> collect

    # Replace in eqns
    out = replace(eqn, r"sin\((\w+)\)" => s"s\1") |>
        x -> replace(x, r"cos\((\w+)\)" => s"c\1")

    # Add to replacements
    for m ∈ [m1; m2]
        lhs = m.captures[1]
        rhs = m.match
        prefix = rhs[1]
        eqn = prefix * lhs * " = " * rhs

        eqn ∉ replacements && push!(replacements, eqn)
    end

    return out
end

for variable ∈ io_all
    # Unwinds Z assignments and simplify
    out = replace_rhs(variable) |> replaceSinCos
    push!(io_raw, out)
end

# Sort outputs
replacements = sort!(replacements) |> x -> reshape(x, :, 2)
sort!(parameters)


## Create files
msg = "# File generated automatically from $fname"
n = length(rhs)
# Convert array assignments to individual variables (for StaticArrays)
function convertArrToVar(x)
    out = similar(x)
    for i ∈ eachindex(x)
        out[i] = replace(x[i], r"(\w+)\[(\d),*(\d)\]" => s"\1\2\3")
        out[i] = replace(out[i], r"(\w+)\[(\d)\]" => s"\1\2")
    end
    
    return out
end
# Find which parameters needed for function
function functionNeeds(funcs, vars)
    out = String[]
    for var ∈ vars
        any(occursin.(var, funcs)) && var ∉ out && push!(out, var)
    end

    return out
end
# Evaluate spline at t from Dierck.jl
evalSpline(p) = "$p = $p(t)"
# Get derivative of of function from Dierckx.jl
getDeriv(p, n) = "$(p * "p"^n) = derivative($p, t, $n)"

# eom function
eomS =
"$msg
# Equations of motion for model
function eom(u, p, t)
    @unpack $(join(parameters, ", ")), $(join(specified_functions, ", ")), z = p 
    @inbounds $(join(variables, ", ")) = u

    $(join(constants[occursin.(r"[q,u]\w+ = ", constants)], "\n\t"))

    $(join([(getDeriv(f, 2)) for f ∈ specified_functions], "\n\t"))
    oxp = derivative(ox, t, 1)
    oyp = derivative(oy, t, 1)
    ozp = derivative(oz, t, 1)
    ox = ox(t)
    oy = oy(t)
    oz = oz(t)

    # Torques
    r1 = RotXYZ(q1 - eqx, q2 - eqy, q3 - eqz); r2 = RotXYZ(q4, q5, q6)
    atorx, atory, atorz = -ka * rotation_angle(r1) * rotation_axis(r1) .- ba * [u1,u2,u3]
    btorx, btory, btorz = -kb * rotation_angle(r2) * rotation_axis(r2) .- bb * [u4,u5,u6]

    $(join(z_qp, "\n\t"))

    $(join(convertArrToVar(coef), "\n\t"))
    $(join(convertArrToVar(rhs), "\n\t"))
    
    coef = @SMatrix [$(begin
        [["coef$j$i" for i ∈ 1:n] |> x -> join(x, ' ') for j ∈ 1:n] |> 
        x -> join(x, "\n\t\t\t\t\t ")
    end)]

    rhs = @SVector [$(["rhs$i" for i ∈ 1:n] |> x -> join(x, ", "))]

    @inbounds $(join(variables[occursin.(r"^u", variables)] .* "p", ", ")) = coef \\ rhs

    @SVector [$(["$(v)p" for v ∈ variables] |> x -> join(x, ", "))]
end  
"
open("model/eom.jl", "w") do io
    write(io, eomS)
end


# creatModel function, allows z array to be visible to eom function
createModelS =
"$( eomS )

# Create model
function createModel(fname)
    
    # Parameters
    px, py, pz, ox, oy, oz,  time, base, mid, tip, initconds, la, lb, tspan = getdata(fname)
	deleteat!(initconds, 10) # Remove u4 
    p = initialise_parameters(px, py, pz, ox, oy, oz, la, lb)

    # Initial conditions
    u₀ = SVector{11,Float64}(initconds)

    prob = ODEProblem(eom, u₀, tspan, p)

    return prob, p, time, base, mid, tip

end
"
open("model/createModel.jl", "w") do io
    write(io, createModelS)
end

## functions function for output quantities
functionsS =
"$msg
"
for f ∈ io_raw
    m = match(r"^(\w+) = (.*)", f)
    name = m.captures[1]
    name == "te" && continue
    eq = m.captures[2]

    f_name = "function $(name)(sol, t)\n"
    # f_param = [functionNeeds(f, parameters); functionNeeds(f, specified_functions)] |> 
    #     x -> "\t@unpack $(join(x, ", ")) = sol.prob.p\n"
    f_param = "\t@unpack ba, bb, eqx, eqy, eqz, g, ixa, ixb, iya, iyb, iza, izb, ka, kb, la, lao, lb, lbo, ma, mb, ox, oy, oz, px, py, pz = sol.prob.p\n"
    f_var = "\t@inbounds $(join(variables, ", ")) = sol(t)\n\n"
    f_cons = "\t$(join(constants[occursin.(r"[q,u]\w+ = ", constants)], "\n\t"))\n"
    f_deriv = "\t$( join(getDeriv.(specified_functions, 1), "\n\t") )\n"
    f_evals = "\t$( join(evalSpline.(specified_functions), "\n\t") )\n\n"
    f_repl = "\t$(join([join(replacements[:,i], "; ") for i ∈ 1:2], "\n\t")) \n\n"
    f_eq = "\treturn $(eq)\n"

    global functionsS *= f_name * f_param * f_var * f_deriv * f_evals * f_repl * f_eq * "end\n\n"

    f2 = "$name(sol) = [$name(sol, t) for t ∈ sol.t]\n\n"
    functionsS *= f2
end

te_fun = "te(sol, t) = ke(sol, t) + pe(sol, t)\n"
te2_fun = "te(sol) = [te(sol, t) for t ∈ sol.t]\n"
functionsS *= te_fun * te2_fun

open("model/functions.jl", "w") do io
    write(io, functionsS)
end


parametersS =
"$msg
# Parameters struct
@with_kw mutable struct Params{T} @deftype T
    z::Vector{Float64} = Vector{Float64}(undef, $(z_array))
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
    $( join(specified_functions, "::Union{Float64,Spline1D}\n\t") )::Union{Float64,Spline1D}
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
    @unpack $( join(functionNeeds(constants, parameters), ", ") ), z = p

    # Evaluate constants
    $( join(constants, "\n\t") )

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
    $( join(constants, "\n\t") )

    return  @pack! p = ixa, ixb, iya, iyb, iza, izb, ma, mb, ka, kb, ba, bb, eqx, eqy, eqz
end
"

open("parameters.jl","w") do io
    write(io, parametersS)
end


# spec = functionNeeds(f, specified)
# f_deriv = ["\t$(getDeriv(s[1:end - 1], 1))" for s ∈ spec if occursin(r"p$", s)] |>
#     x -> join(x, "\n")
# isempty(f_deriv) || (f_deriv *= "\n\n")
# evals = functionNeeds(f, specified_functions)
# f_evals = ["\t$(evalSpline(e))" for e ∈ evals] |> x -> join(x, '\n')
# isempty(f_evals) || (f_evals *= "\n\n")
# f_repl = "\t$(join([join(replacements[:,i], "; ") for i ∈ 1:2], "\n\t")) \n\n"

