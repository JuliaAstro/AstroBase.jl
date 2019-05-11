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

        write(f, """
const VSOP_$(name) = ((($(name)_X_0_AMP, $(name)_X_0_PHS, $(name)_X_0_FRQ),
                       ($(name)_X_1_AMP, $(name)_X_1_PHS, $(name)_X_1_FRQ),
                       ($(name)_X_2_AMP, $(name)_X_2_PHS, $(name)_X_2_FRQ),
                       ($(name)_X_3_AMP, $(name)_X_3_PHS, $(name)_X_3_FRQ),
                       ($(name)_X_4_AMP, $(name)_X_4_PHS, $(name)_X_4_FRQ),
                       ($(name)_X_5_AMP, $(name)_X_5_PHS, $(name)_X_5_FRQ)),
                      (($(name)_Y_0_AMP, $(name)_Y_0_PHS, $(name)_Y_0_FRQ),
                       ($(name)_Y_1_AMP, $(name)_Y_1_PHS, $(name)_Y_1_FRQ),
                       ($(name)_Y_2_AMP, $(name)_Y_2_PHS, $(name)_Y_2_FRQ),
                       ($(name)_Y_3_AMP, $(name)_Y_3_PHS, $(name)_Y_3_FRQ),
                       ($(name)_Y_4_AMP, $(name)_Y_4_PHS, $(name)_Y_4_FRQ),
                       ($(name)_Y_5_AMP, $(name)_Y_5_PHS, $(name)_Y_5_FRQ)),
                      (($(name)_Z_0_AMP, $(name)_Z_0_PHS, $(name)_Z_0_FRQ),
                       ($(name)_Z_1_AMP, $(name)_Z_1_PHS, $(name)_Z_1_FRQ),
                       ($(name)_Z_2_AMP, $(name)_Z_2_PHS, $(name)_Z_2_FRQ),
                       ($(name)_Z_3_AMP, $(name)_Z_3_PHS, $(name)_Z_3_FRQ),
                       ($(name)_Z_4_AMP, $(name)_Z_4_PHS, $(name)_Z_4_FRQ),
                       ($(name)_Z_5_AMP, $(name)_Z_5_PHS, $(name)_Z_5_FRQ)))

const VSOP_$(name)_NUM = ($(name)_X_NUM, $(name)_Y_NUM, $(name)_Z_NUM)
""")

    end
    empty!(coeffs[1])
    empty!(coeffs[2])
    empty!(coeffs[3])
end
