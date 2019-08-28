using SPICE:
    SpiceError,
    bodfnd,
    bodn2c,
    bodvcd,
    eul2m,
    eul2xf,
    furnsh,
    tipbod,
    tisbod

furnsh(joinpath(@__DIR__, "..", "gen", "gm_de431.tpc"),
       joinpath(@__DIR__, "..", "gen", "pck00010.tpc"))

@testset "Celestial Bodies" begin
    @testset for body_name in AstroBase.Bodies.ALL_NAMES
        body = @eval $(Symbol(body_name))()
        id = naifid(body)
        name = string(body)
        name = name == "Luna" ? "Moon" : name

        @testset "NAIF ID" begin
            @test id == bodn2c(name)
            @test from_naifid(naifid(body)) == body
        end

        if bodfnd(id, "GM")
            @testset "Gravitational parameter" begin
                μ = grav_param(body)
                @test μ ≈ bodvcd(id, "GM")[1]
            end
        end

        if bodfnd(id, "PM")
            @testset "Rotational elements" begin
                ep = TDBEpoch(2019, 5, 6)
                et = value(seconds(j2000(ep)))
                eulang = [reverse(euler_angles(body, ep))...; reverse(euler_rates(body, ep))...]
                p_act = eul2m(eulang[1:3]..., 3, 1, 3)
                p_exp = tipbod("J2000", id, et)
                @testset for i in eachindex(p_act, p_exp)
                    @test p_act[i] ≈ p_exp[i] atol=1e-8
                end
                s_act = eul2xf(eulang, 3, 1, 3)
                s_exp = tisbod("J2000", id, et)
                @testset for i in eachindex(s_act, s_exp)
                    @test s_act[i] ≈ s_exp[i] atol=1e-8
                end
            end
        end
    end
    @test AstroBase.Bodies.path_ids(luna, earth) == [301, 3, 399]
    @test AstroBase.Bodies.path_ids(luna, io) == [301, 3, 0, 5, 501]
end
