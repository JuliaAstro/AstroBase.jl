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

        nut80_exp = ERFA.nut80(jd...)
        nut80_act = nutation(iau1980, ep)
        @test nut80_act[1] ≈ nut80_exp[1]
        @test nut80_act[2] ≈ nut80_exp[2]

        nut00a_exp = ERFA.nut00a(jd...)
        nut00a_act = nutation(iau2000a, ep)
        @test nut00a_act[1] ≈ nut00a_exp[1]
        @test nut00a_act[2] ≈ nut00a_exp[2]

        nut00b_exp = ERFA.nut00b(jd...)
        nut00b_act = nutation(iau2000b, ep)
        @test nut00b_act[1] ≈ nut00b_exp[1]
        @test nut00b_act[2] ≈ nut00b_exp[2]

        nut06_exp = ERFA.nut06a(jd...)
        nut06_act = nutation(iau2006, ep)
        @test nut06_act[1] ≈ nut06_exp[1]
        @test nut06_act[2] ≈ nut06_exp[2]
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
    @test numat(0.7, 1.4, 1.3) ≈ ERFA.numat(0.7, 1.4, 1.3)

    let (ap, ao) = precession_rate_part_of_nutation(2.4578265e6, 0.30434616919175345)
        (ep, eo) = ERFA.pr00(2.4578265e6, 0.30434616919175345)
        @test ap ≈ ep
        @test ao ≈ eo
    end

    @test greenwich_mean_sidereal_time00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test greenwich_mean_sidereal_time06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test greenwich_mean_sidereal_time82(2.4578265e6, 0.30434616919175345) ≈ ERFA.gmst82(2.4578265e6, 0.30434616919175345)
    @test greenwich_mean_sidereal_time82(0.30434616919175345, 2.4578265e6) ≈ ERFA.gmst82(0.30434616919175345, 2.4578265e6)
    let (b1, p1, bp1) = bias_precession_matrix_00(2.4578265e6, 0.30434616919175345)
        (b2, p2, bp2) = ERFA.bp00(2.4578265e6, 0.30434616919175345)
        @test b1 ≈ b2
        @test p1 ≈ p2
        @test bp1 ≈ bp2
    end

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
    @test nutation_matrix80(2.4578265e6, 0.30434616919175345) ≈ ERFA.nutm80(2.4578265e6, 0.30434616919175345)

    let (a1, b1, c1, d1, e1, f1) = precession_nutation00(2.4578265e6, 0.30434616919175345, 0.2, 0.2)
        (a2, b2, c2, d2, e2, f2) = ERFA.pn00(2.4578265e6, 0.30434616919175345, 0.2, 0.2)
        @test a1 ≈ a2
        @test b1 ≈ b2
        @test c1 ≈ c2
        @test d1 ≈ d2
        @test e1 ≈ e2
        @test f1 ≈ f2
    end

    let (a1, b1, c1) = precession_nutation_a00(2.4578265e6, 0.30434616919175345)
        (a2, b2, c2, d2, e2, f2, g2, h2) = ERFA.pn00a(2.4578265e6, 0.30434616919175345)
        @test a1 ≈ a2
        @test b1 ≈ b2
        @test c1[1] ≈ c2
        @test c1[2] ≈ d2
        @test c1[3] ≈ e2
        @test c1[4] ≈ f2
        @test c1[5] ≈ g2
        @test c1[6] ≈ h2
    end

    let (a1, b1, c1) = precession_nutation_b00(2.4578265e6, 0.30434616919175345)
        (a2, b2, c2, d2, e2, f2, g2, h2) = ERFA.pn00b(2.4578265e6, 0.30434616919175345)
        @test a1 ≈ a2
        @test b1 ≈ b2
        @test c1[1] ≈ c2
        @test c1[2] ≈ d2
        @test c1[3] ≈ e2
        @test c1[4] ≈ f2
        @test c1[5] ≈ g2
        @test c1[6] ≈ h2
    end

    @test precession_nutation_matrix_a00(2.4578265e6, 0.30434616919175345) ≈ ERFA.pnm00a(2.4578265e6, 0.30434616919175345)
    @test precession_nutation_matrix_b00(2.4578265e6, 0.30434616919175345) ≈ ERFA.pnm00b(2.4578265e6, 0.30434616919175345)
    @test precession_nutation_matrix_a06(2.4578265e6, 0.30434616919175345) ≈ ERFA.pnm06a(2.4578265e6, 0.30434616919175345)

    @test nutation_matrix_day_a00(2.4578265e6, 0.30434616919175345) ≈ ERFA.num00a(2.4578265e6, 0.30434616919175345)
    @test nutation_matrix_day_b00(2.4578265e6, 0.30434616919175345) ≈ ERFA.num00b(2.4578265e6, 0.30434616919175345)

    @test nutation_matrix_day(2.4578265e6, 0.30434616919175345) ≈ ERFA.num06a(2.4578265e6, 0.30434616919175345)
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

