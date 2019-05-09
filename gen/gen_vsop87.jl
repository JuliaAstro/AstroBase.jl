files = Dict("earth" => "VSOP87E.ear",
             "jupiter" => "VSOP87E.jup",
             "mars" => "VSOP87E.mar",
             "mercury" => "VSOP87E.mer",
             "neptune" => "VSOP87E.nep",
             "saturn" => "VSOP87E.sat",
             "sun" => "VSOP87E.sun",
             "uranus" => "VSOP87E.ura",
             "venus" => "VSOP87E.ven")

coeffs = (x=[], y=[], z=[])

for (name, f) in files
    lines = open(readlines, f)
    for line in lines
        startswith(line, " VSOP") && continue

        idx = parse(Int, line[4:4])
        order = parse(Int, line[5:5]) + 1
        if length(coeffs[idx]) < order
            push!(coeffs[idx], (amplitude=[], phase=[], frequency=[]))
        end
        push!(coeffs[idx][order][:amplitude], parse(Float64, line[80:97]))
        push!(coeffs[idx][order][:phase], parse(Float64, line[98:111]))
        push!(coeffs[idx][order][:frequency], parse(Float64, line[112:131]))
    end

    open("vsop87_$name.jl", "w") do f
        for (idx, coord) in enumerate(["x", "y", "z"])
            write(f, "const $(name)_$(coord) = (\n")
            for (order, co) in enumerate(coeffs[idx])
                write(f, "$(repeat(" ", 4))(\n$(repeat(" ", 8))[\n")
                for c in co[1]
                    write(f, "$(repeat(" ", 12))$c,\n")
                end
                write(f, "$(repeat(" ", 8))],\n$(repeat(" ", 8))[\n")
                for c in co[2]
                    write(f, "$(repeat(" ", 12))$c,\n")
                end
                write(f, "$(repeat(" ", 8))],\n$(repeat(" ", 8))[\n")
                for c in co[3]
                    write(f, "$(repeat(" ", 12))$c,\n")
                end
                write(f, "$(repeat(" ", 8))]\n$(repeat(" ", 4))),\n")
            end
            write(f, ")\n\n")
        end
    end
    empty!(coeffs[1])
    empty!(coeffs[2])
    empty!(coeffs[3])
end
