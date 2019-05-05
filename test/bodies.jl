using SPICE

furnsh(joinpath(@__DIR__, "..", "gen", "gm_de431.tpc"),
       joinpath(@__DIR__, "..", "gen", "pck00010.tpc"))

@testset "Celestial Bodies" begin
    for body in AstroBase.Bodies.ALL_NAMES
        typ = @eval $(Symbol(body))
        id = naifid(typ())
        name = string(typ())
        name = name == "Luna" ? "Moon" : name
        @test id == bodn2c(name)
        @test from_naifid(naifid(typ())) == typ()
        try
            μ = grav_param(typ())
            @test μ == bodvcd(id, "GM")[1]
        catch err
            err isa MethodError || rethrow(err)
        end
    end
    @test AstroBase.Bodies.path_ids(luna, earth) == [301, 3, 399]
    @test AstroBase.Bodies.path_ids(luna, io) == [301, 3, 0, 5, 501]
end
