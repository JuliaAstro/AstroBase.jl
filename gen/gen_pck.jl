using SPICE
using Statistics: mean

file = "pck00010.tpc"
furnsh(file)

ids = Set([])
re = r"BODY(?<id>[0-9]+)_(RADII|PM)\s+="
lines = open(readlines, file)
for line in lines
    s = string(line)
    if occursin(re, s)
        m = match(re, s)
        push!(ids, m["id"])
    end
end

template(body, func, val) = """
@inline $func(::Type{Float64}, ::$body) = $val
@inline $func(::$body) = $func(Float64, $body())
"""

open("pck.jl", "w") do f
    write(f, "# Automatically generated by $(basename(@__FILE__)), do not edit!\n\n")
    for id in ids
        body = replace(titlecase(bodc2n(parse(Int, id))), r"\s"=>"")
        body = body == "Moon" ? "Luna" : body

        # Radii
        try
            r = bodvrd(id, "RADII")
            subplan, along, polar = r
            write(f, template(body, "subplanetary_radius", subplan))
            write(f, template(body, "along_orbit_radius", along))
            if subplan == along
                write(f, template(body, "equatorial_radius", subplan))
            end
            write(f, template(body, "polar_radius", polar))
            write(f, template(body, "mean_radius", mean(r)))
        catch err
            err isa SpiceError || rethrow(err)
        end

        functions = ("right_ascension", "declination", "rotation")
        names = ("RA", "DEC", "PM")

        try
            for (func, name) in zip(functions, names)
                val = zeros(3)
                sval = deg2rad.(bodvrd(id, name == "PM" ? name : "POLE_$name"))
                val[1:length(sval)] .+= sval
                val_np = Float64[]
                try
                    n = length(bodvrd(id[1:1], "NUT_PREC_ANGLES")) ÷ 2
                    tmp = zeros(n)
                    sval = deg2rad.(bodvrd(id, "NUT_PREC_$name"))
                    tmp[1:length(sval)] .+= sval
                    append!(val_np, tmp)
                catch err
                    err isa SpiceError || rethrow(err)
                end
                data = (val..., Tuple(val_np))
                write(f, template(body, "$(func)_coeffs", data))
            end

            # Nutation / Precession
            if length(id) != 3
                write(f, template(body, "nutation_precession_coeffs", ((),())))
                continue
            end
            try
                np = deg2rad.(reshape(bodvrd(id[1:1], "NUT_PREC_ANGLES"), 2, :))
                theta0 = Tuple(np[1,:])
                theta1 = Tuple(np[2,:])
                write(f, template(body, "nutation_precession_coeffs", (theta0, theta1)))
            catch err
                err isa SpiceError || rethrow(err)
                write(f, template(body, "nutation_precession_coeffs", ((),())))
            end
        catch err
            err isa SpiceError || rethrow(err)
        end
    end
end
