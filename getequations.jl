# Reads fortran file output from Autolev and rewrites equations in EQNS and IO in Julia

# Read whole file as string
lines = read("model/dogtail.f", String)

# Qps'
msqp = eachmatch(r"Q[1-9][0-9]?p =\s+.*(\n\s+&.*){0,10}", lines)
qp = []
for mqp ∈ msqp
    # Get full assignment, remove indexing
    out = mqp.match

    # Remove continuation characters
    out = replace(out, r"(\r\n\s+&|\r)" => s"")

    # Replace exponentials
    out = replace(out, "**" => "^")

    # Make lowercase
    out = lowercase(out)

    # Replace sin/cos with s/c
    out = replace(out, r"sin\(\w(\d)\)" => s"s\1")
    out = replace(out, r"cos\(\w(\d)\)" => s"c\1")

    # Remove d0
    out = replace(out, "d0" => "")

    # println(out)
    push!(qp, out)

end

# COEF
mscoef = eachmatch(r"COEF\((\d),(\d)\)\s+=((.*\n\s+&){0,10}(.*))", lines)
coef = []
for mcoef ∈ mscoef
    # Get full assignment, remove indexing
    out = "COEF$(mcoef.captures[1])$(mcoef.captures[2]) = $(mcoef.captures[3])"

    # Remove continuation characters
    out = replace(out, r"(\r\n\s+&|\r)" => s"")

    # Replace exponentials
    out = replace(out, "**" => "^")

    # Make lowercase
    out = lowercase(out)

    # Replace sin/cos with s/c
    out = replace(out, r"sin\(\w(\d)\)" => s"s\1")
    out = replace(out, r"cos\(\w(\d)\)" => s"c\1")

    # Remove d0
    out = replace(out, "d0" => "")

    # println(out)
    push!(coef, out)

end

# RHS
msrhs = eachmatch(r"RHS\((\d)\)\s+=((.*\n\s+&){0,200}(.*))", lines)
rhs = []
for mrhs ∈ msrhs
    # Get full assignment, remove indexing
    out = "RHS$(mrhs.captures[1]) = $(mrhs.captures[2])"
    
    # Remove continuation characters
    out = replace(out, r"(\r\n\s+&|\r)" => s"")

    # Replace exponentials
    out = replace(out, "**" => "^")

    # Make lowercase
    out = lowercase(out)

    # Replace sin/cos with s/c
    out = replace(out, r"sin\(\w(\d)\)" => s"s\1")
    out = replace(out, r"cos\(\w(\d)\)" => s"c\1")

    # Remove d0
    out = replace(out, "d0" => "")

    # println(out)
    push!(rhs, out)
    
end

# Points
mspoints = eachmatch(r"P[2,3][X-Z] =\s+.*(\n\s+&.*){0,10}", lines)
points = []
for mpoint ∈ mspoints
    # Get full string
    out = mpoint.match

    # Remove continuation characters
    out = replace(out, r"(\r\n\s+&|\r)" => s"")

    # Replace exponentials
    out = replace(out, "**" => "^")

    # Make lowercase
    out = lowercase(out)

    # Replace sin/cos with s/c
    out = replace(out, r"sin\(\w(\d)\)" => s"s\1")
    out = replace(out, r"cos\(\w(\d)\)" => s"c\1")

    # Remove d0
    out = replace(out, "d0" => "")
    
    # println(out)
    push!(points, out)
end

# Energy
msenergy = eachmatch(r"[K,P]E =\s+.*(\n\s+&.*){0,100}", lines)
energy = []
for menergy ∈ msenergy
    # Get full string
    out = menergy.match

    # Remove continuation characters
    out = replace(out, r"(\r\n\s+&|\r)" => s"")

    # Replace exponentials
    out = replace(out, "**" => "^")

    # Make lowercase
    out = lowercase(out)

    # Replace sin/cos with s/c
    # out = replace(out, r"sin\(\w(\d)\)" => s"s\1")
    # out = replace(out, r"cos\(\w(\d)\)" => s"c\1")

    # Remove d0
    out = replace(out, "d0" => "")

    # println(out)
    push!(energy, out)
end

# Momenta
msmomenta = eachmatch(r"AMOM[X-Z] =\s+.*(\n\s+&.*){0,100}", lines)
momentum =  []
for mmomenta ∈ msmomenta
    # Get full string
    out = mmomenta.match

    # Remove continuation characters
    out = replace(out, r"(\r\n\s+&|\r)" => s"")

    # Replace exponentials
    out = replace(out, "**" => "^")

    # Make lowercase
    out = lowercase(out)

    # Replace sin/cos with s/c
    # out = replace(out, r"sin\(\w(\d)\)" => s"s\1")
    # out = replace(out, r"cos\(\w(\d)\)" => s"c\1")

    # Remove d0
    out = replace(out, "d0" => "")

    # println(out)
    push!(momentum, out)
end


# Write to output file
open("model/equations.jl", "w") do io
    [println(io, eq) for eq ∈ qp]
    println(io, "")
    [println(io, eq) for eq ∈ coef]
    [println(io, eq) for eq ∈ rhs]
    println(io, "")
    [println(io, eq) for eq ∈ points]
    println(io, "")
    [println(io, eq) for eq ∈ energy]
    println(io, "")
    [println(io, eq) for eq ∈ momentum]

end