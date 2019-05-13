files = Dict("EARTH" => "VSOP87E.ear",
             "JUPITER" => "VSOP87E.jup",
             "MARS" => "VSOP87E.mar",
             "MERCURY" => "VSOP87E.mer",
             "NEPTUNE" => "VSOP87E.nep",
             "SATURN" => "VSOP87E.sat",
             "SUN" => "VSOP87E.sun",
             "URANUS" => "VSOP87E.ura",
             "VENUS" => "VSOP87E.ven")

coeffs = (x=[], y=[], z=[])

indent = repeat(" ", 4)

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
        for (idx, coord) in enumerate(["X", "Y", "Z"])
            write(f, "const $(name)_$(coord)_NUM = (\n")
            for (order, co) in enumerate(coeffs[idx])
                write(f, "$indent$(length(co[1])),\n")
            end
            write(f, ")\n\n")
            for (i, co) in enumerate(coeffs[idx])
                order = i - 1
                write(f, "const $(name)_$(coord)_$(order)_AMP = [\n")
                for c in co[1]
                    write(f, "$indent$c,\n")
                end
                write(f, "]\n\n")
                write(f, "const $(name)_$(coord)_$(order)_PHS = [\n")
                for c in co[2]
                    write(f, "$indent$c,\n")
                end
                write(f, "]\n\n")
                write(f, "const $(name)_$(coord)_$(order)_FRQ = [\n")
                for c in co[3]
                    write(f, "$indent$c,\n")
                end
                write(f, "]\n\n")
            end
        end

        write(f, "const VSOP_$(name) = ((\n")
        for (idx, coord) in enumerate(["X", "Y", "Z"])
            for i in 0:length(coeffs[idx]) - 1
                write(f, "$indent($(name)_$(coord)_$(i)_AMP, $(name)_$(coord)_$(i)_PHS, $(name)_$(coord)_$(i)_FRQ),\n")
            end
            if idx in (1, 2)
                write(f, "), (\n")
            end
        end
        write(f, "))\n\n")
        write(f, "const VSOP_$(name)_NUM = ($(name)_X_NUM, $(name)_Y_NUM, $(name)_Z_NUM)\n\n")
    end
    empty!(coeffs[1])
    empty!(coeffs[2])
    empty!(coeffs[3])
end
