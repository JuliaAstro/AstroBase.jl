@testset "Celestial Bodies" begin
    let bodies = ["Sun"]
        append!(bodies, AstroBase.Bodies.PLANET_NAMES)
        append!(bodies, ["Luna", "Phobos", "Deimos"])
        append!(bodies, AstroBase.Bodies.JUPITER_SATELLITE_NAMES)
        append!(bodies, AstroBase.Bodies.SATURN_SATELLITE_NAMES)
        append!(bodies, AstroBase.Bodies.URANUS_SATELLITE_NAMES)
        append!(bodies, AstroBase.Bodies.NEPTUNE_SATELLITE_NAMES)
        append!(bodies, AstroBase.Bodies.MINOR_BODY_NAMES)
        push!(bodies, "Pluto")
        append!(bodies, AstroBase.Bodies.PLUTO_SATELLITE_NAMES)
        for body in bodies
            typ = @eval Symbol($body)
            @test string(typ) == body
        end
    end
end
