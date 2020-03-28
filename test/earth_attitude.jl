using AstroBase
using Test

import ERFA

@testset "Earth Attitude" begin
    @testset "Obliquity" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))

        @test obliquity(iau1980, ep) ≈ ERFA.obl80(jd...)
        @test obliquity(iau2006, ep) ≈ ERFA.obl06(jd...)
    end
    @testset "Nutation" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))

        @testset "nut80" begin
            nut80_exp = ERFA.nut80(jd...)
            nut80_act = nutation(iau1980, ep)
            @test nut80_act[1] ≈ nut80_exp[1]
            @test nut80_act[2] ≈ nut80_exp[2]
        end
        @testset "nut00a" begin
            nut00a_exp = ERFA.nut00a(jd...)
            nut00a_act = nutation(iau2000a, ep)
            @test nut00a_act[1] ≈ nut00a_exp[1]
            @test nut00a_act[2] ≈ nut00a_exp[2]
        end
        @testset "nut00b" begin
            nut00b_exp = ERFA.nut00b(jd...)
            nut00b_act = nutation(iau2000b, ep)
            @test nut00b_act[1] ≈ nut00b_exp[1]
            @test nut00b_act[2] ≈ nut00b_exp[2]
        end
        @testset "nut06a" begin
            nut06a_exp = ERFA.nut06a(jd...)
            nut06a_act = nutation(iau2006a, ep)
            @test nut06a_act[1] ≈ nut06a_exp[1]
            @test nut06a_act[2] ≈ nut06a_exp[2]
        end
        @testset "numat" begin
            ϵ = obliquity(iau2006, ep)
            δψ, δϵ = nutation(iau2006a, ep)
            numat_exp = ERFA.numat(ϵ, δψ, δϵ)
            numat_act = nutation_matrix(ϵ, δψ, δϵ)
            @testset for i in eachindex(numat_act, numat_exp)
                @test numat_act[i] ≈ numat_exp[i]
            end
        end
        @testset "nutm80" begin
            nutm80_exp = ERFA.nutm80(jd...)
            nutm80_act = nutation_matrix(iau1980, ep)
            @testset for i in eachindex(nutm80_act, nutm80_exp)
                @test nutm80_act[i] ≈ nutm80_exp[i]
            end
        end
        @testset "num00a" begin
            num00a_exp = ERFA.num00a(jd...)
            num00a_act = nutation_matrix(iau2000a, ep)
            @testset for i in eachindex(num00a_act, num00a_exp)
                @test num00a_act[i] ≈ num00a_exp[i]
            end
        end
        @testset "num00b" begin
            num00b_exp = ERFA.num00b(jd...)
            num00b_act = nutation_matrix(iau2000b, ep)
            @testset for i in eachindex(num00b_act, num00b_exp)
                @test num00b_act[i] ≈ num00b_exp[i]
            end
        end
        @testset "num06a" begin
            num06a_exp = ERFA.num06a(jd...)
            num06a_act = nutation_matrix(iau2006a, ep)
            @testset for i in eachindex(num06a_act, num06a_exp)
                @test num06a_act[i] ≈ num06a_exp[i]
            end
        end
    end
    @testset "Precession" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))

        @testset "bi00" begin
            bi00_exp = ERFA.bi00()
            bi00_act = bias(iau2000)
            @testset for i in eachindex(bi00_act, bi00_exp)
                @test bi00_act[i] ≈ bi00_exp[i]
            end
        end
        @testset "pr00" begin
            pr00_exp = ERFA.pr00(jd...)
            pr00_act = precession(iau2000, ep)
            @test pr00_act[1] ≈ pr00_exp[1]
            @test pr00_act[2] ≈ pr00_exp[2]

            bp00_exp = ERFA.bp00(jd...)
            bp00_act = bias_precession_matrix(iau2000, ep)
            @testset for i in eachindex(bp00_act)
                act = bp00_act[i]
                exp = bp00_exp[i]
                @testset for j in eachindex(act, exp)
                    @test act[j] ≈ exp[j]
                end
            end
        end
    end
    @testset "Precession-Nutation" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))

        @testset "pn00" begin
            nut = nutation(iau2000a, ep)
            pn00_exp = ERFA.pn00(jd..., nut...)
            pn00_act = precession_nutation(iau2000, ep, nut...)
            @testset for i in eachindex(pn00_act, pn00_exp)
                if i == 1
                    @test pn00_act[i] ≈ pn00_exp[i]
                else
                    act = pn00_act[i]
                    exp = pn00_exp[i]
                    @testset for j in eachindex(act, exp)
                        @test act[j] ≈ exp[j]
                    end
                end
            end
        end
        @testset "pn00a" begin
            pn00a_exp = ERFA.pn00a(jd...)
            pn00a_act = precession_nutation(iau2000a, ep)
            @testset for i in eachindex(pn00a_act, pn00a_exp)
                if i in 1:3
                    @test pn00a_act[i] ≈ pn00a_exp[i]
                else
                    act = pn00a_act[i]
                    exp = pn00a_exp[i]
                    @testset for j in eachindex(act, exp)
                        @test act[j] ≈ exp[j]
                    end
                end
            end
        end
        @testset "pn00b" begin
            pn00b_exp = ERFA.pn00b(jd...)
            pn00b_act = precession_nutation(iau2000b, ep)
            @testset for i in eachindex(pn00b_act, pn00b_exp)
                if i in 1:3
                    @test pn00b_act[i] ≈ pn00b_exp[i]
                else
                    act = pn00b_act[i]
                    exp = pn00b_exp[i]
                    @testset for j in eachindex(act, exp)
                        @test act[j] ≈ exp[j]
                    end
                end
            end
        end
        @testset "pnm00a" begin
            pnm00a_exp = ERFA.pnm00a(jd...)
            pnm00a_act = precession_nutation_matrix(iau2000a, ep)
            @testset for i in eachindex(pnm00a_act, pnm00a_exp)
                @test pnm00a_act[i] ≈ pnm00a_exp[i]
            end
        end
        @testset "pnm00b" begin
            pnm00b_exp = ERFA.pnm00b(jd...)
            pnm00b_act = precession_nutation_matrix(iau2000b, ep)
            @testset for i in eachindex(pnm00b_act, pnm00b_exp)
                @test pnm00b_act[i] ≈ pnm00b_exp[i]
            end
        end
    end
    @test earth_rotation_angle(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)

    let (w1, x1, y1, z1) = precession_fukushima_williams06(2.4578265e6, 0.30434616919175345)
        (w2, x2, y2, z2) = ERFA.pfw06(2.4578265e6, 0.30434616919175345)
        @test w1 ≈ w2
        @test x1 ≈ x2
        @test y1 ≈ y2
        @test z1 ≈ z2
    end

    # Celestial to intermediate frame of date
    @test celestial_to_intermediate(0.2, 0.2, 0.2) ≈ ERFA.c2ixys(0.2, 0.2, 0.2)

    # Polar motion
    @test polar_motion(30, 30, 30) ≈ ERFA.pom00(30, 30, 30)
    @test polar_motion(20, 30, 50) ≈ ERFA.pom00(20, 30, 50)

    @test tio_locator(2.4578265e6, 0.30434616919175345) ≈ ERFA.sp00(2.4578265e6, 0.30434616919175345)

    @test AstroBase.EarthAttitude.mean_anomaly(luna, 1)  ≈ ERFA.fal03(1)
    @test AstroBase.EarthAttitude.mean_anomaly(sun, 1) ≈ ERFA.falp03(1)
    @test mean_longitude_minus_lan(luna, 1) ≈ ERFA.faf03(1)
    @test mean_elongation(luna, 1) ≈ ERFA.fad03(1)
    @test mean_longitude_ascending_node(luna, 1) ≈ ERFA.faom03(1)
    @test mean_longitude(mercury, 1) ≈ ERFA.fame03(1)
    @test mean_longitude(venus, 1) ≈ ERFA.fave03(1)
    @test mean_longitude(earth, 1) ≈ ERFA.fae03(1)
    @test mean_longitude(mars, 1) ≈ ERFA.fama03(1)
    @test mean_longitude(jupiter, 1) ≈ ERFA.faju03(1)
    @test mean_longitude(saturn, 1) ≈ ERFA.fasa03(1)
    @test mean_longitude(uranus, 1) ≈ ERFA.faur03(1)
    @test mean_longitude(neptune, 1) ≈ ERFA.fane03(1)
    @test general_precession_in_longitude(1) ≈ ERFA.fapa03(1)
    let (x1, y1) = xy06(2.4578265e6, 0.30440190993249416)
        (x2, y2) = ERFA.xy06(2.4578265e6, 0.30440190993249416)
        @test x1 ≈ x2
        @test y1 ≈ y2
    end

    @test fukushima_williams_matrix(0.2,0.3,0.5,0.6) ≈ ERFA.fw2m(0.2,0.3,0.5,0.6)

    @test greenwich_mean_sidereal_time00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test greenwich_mean_sidereal_time06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test greenwich_mean_sidereal_time82(2.4578265e6, 0.30434616919175345) ≈ ERFA.gmst82(2.4578265e6, 0.30434616919175345)
    @test greenwich_mean_sidereal_time82(0.30434616919175345, 2.4578265e6) ≈ ERFA.gmst82(0.30434616919175345, 2.4578265e6)

    let a = [0.9999989440476103608 -0.1332881761240011518e-2 -0.5790767434730085097e-3;
         0.1332858254308954453e-2 0.9999991109044505944 -0.4097782710401555759e-4;
         0.5791308472168153320e-3 0.4020595661593994396e-4 0.9999998314954572365]
        @test equation_of_origins(a, 0.3) ≈ ERFA.eors(a, 0.3)
    end

    @test equation_of_equinoxes_complementary_terms(2.4578265e6, 0.30434616919175345) ≈ ERFA.eect00(2.4578265e6, 0.30434616919175345)
    @test equation_of_equinoxes_00(2.4578265e6, 0.30440190993249416, 1.5, 1.7) ≈ ERFA.ee00(2.4578265e6, 0.30440190993249416, 1.5, 1.7)

    @test s00(2.4578265e6, 0.30434616919175345, 0, 0) ≈ ERFA.s00(2.4578265e6, 0.30434616919175345, 0, 0)
    @test s06(2.4578265e6, 0.30434616919175345, 0, 0) ≈ ERFA.s06(2.4578265e6, 0.30434616919175345, 0, 0)

    let t = [0.9999989440476103608 -0.1332881761240011518e-2 -0.5790767434730085097e-3;
            0.1332858254308954453e-2 0.9999991109044505944 -0.4097782710401555759e-4;
            0.5791308472168153320e-3 0.4020595661593994396e-4 0.9999998314954572365]
        @test greenwich_apparent_sidereal_time06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, t) ≈ ERFA.gst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, t)
    end
    @test equation_of_equinoxes_94(2.4578265e6, 0.30434616919175345) ≈ equation_of_equinoxes_94(2.4578265e6, 0.30434616919175345)
end

@testset "Dependencies" begin
    @test s06a(2.4578265e6, 0.30434616919175345) ≈ ERFA.s06a(2.4578265e6, 0.30434616919175345)

    @test equation_of_equinoxes_a00(2.4578265e6, 0.30434616919175345) ≈ ERFA.ee00a(2.4578265e6, 0.30434616919175345)
    @test equation_of_equinoxes_b00(2.4578265e6, 0.30434616919175345) ≈ ERFA.ee00b(2.4578265e6, 0.30434616919175345)

    @test greenwich_mean_sidereal_time_a06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈ ERFA.gst06a(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test equation_of_equinoxes_a06(2.4578265e6, 0.30434616919175345) ≈ ERFA.ee06a(2.4578265e6, 0.30434616919175345)
    @test greenwich_apparent_sidereal_time94(2.4578265e6, 0.30434616919175345) ≈ ERFA.gst94(2.4578265e6, 0.30434616919175345)
    @test greenwich_apparent_sidereal_time_a00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈ ERFA.gst00a(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test greenwich_apparent_sidereal_time_b00(2.4579405e6, 0.0) ≈ ERFA.gst00b(2.4579405e6, 0.0)
    @test greenwich_apparent_sidereal_time_a06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈ ERFA.gst06a(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test celestial_to_intermediate_frame_of_date(2.4578265e6, 0.30434616919175345, 0.2, 0.3) ≈ ERFA.c2ixy(2.4578265e6, 0.30434616919175345, 0.2, 0.3)
    let t = [9.999962358680738e-1 -2.516417057665452e-3 -1.093569785342370e-3;
            2.516462370370876e-3 9.999968329010883e-1 4.006159587358310e-5;
            1.093465510215479e-3 -4.281337229063151e-5 9.999994012499173e-1]
        @test celestial_to_intermediate_matrix(2.4578265e6, 0.30434616919175345, t) ≈ ERFA.c2ibpn(2.4578265e6, 0.30434616919175345, t)
    end
    @test celestial_to_intermediate_matrix_a00(2.4578265e6, 0.30434616919175345) ≈ ERFA.c2i00a(2.4578265e6, 0.30434616919175345)
    @test celestial_to_intermediate_matrix_b00(2.4578265e6, 0.30434616919175345) ≈ ERFA.c2i00b(2.4578265e6, 0.30434616919175345)
    @test celestial_to_terrestrial_matrix(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3, 0.5, 0.7) ≈ ERFA.c2txy(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3, 0.5, 0.7)
    @test celestial_to_terrestrial_a00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3) ≈ ERFA.c2t00a(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3)
    @test celestial_to_terrestrial_b00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3) ≈ ERFA.c2t00b(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3)
    @test c2tpe(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3, 0.5, 0.7) ≈ ERFA.c2tpe(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3, 0.5, 0.7)

    @test s00a(2.4578265e6, 0.30434616919175345) ≈ ERFA.s00a(2.4578265e6, 0.30434616919175345)
    @test s00b(2.4578265e6, 0.30434616919175345) ≈ ERFA.s00b(2.4578265e6, 0.30434616919175345)

    let (x1, y1, z1) = xys00a(2.4578265e6, 0.30434616919175345)
        (x2, y2, z2) = ERFA.xys00a(2.4578265e6, 0.30434616919175345)
        @test x1 ≈ x2
        @test y1 ≈ y2
        @test z1 ≈ z2
    end
    let (x1, y1, z1) = xys00b(2.4578265e6, 0.30434616919175345)
        (x2, y2, z2) = ERFA.xys00b(2.4578265e6, 0.30434616919175345)
        @test x1 ≈ x2
        @test y1 ≈ y2
        @test z1 ≈ z2
    end
end

