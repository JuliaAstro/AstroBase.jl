using AstroBase
using Test

import ERFA

@testset "Earth Attitude" begin
    @testset "Fundamental Arguments" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        t = julian_period(Float64, ep; scale=TDB, unit=centuries)
        longitude = AstroBase.EarthAttitude.Longitude()
        elongation = AstroBase.EarthAttitude.Elongation()
        node = AstroBase.EarthAttitude.AscendingNode()
        atol = 1e-12
        @test fundamental(luna, t) ≈ ERFA.fal03(t) atol=atol
        @test fundamental(sun, t) ≈ ERFA.falp03(t) atol=atol
        @test fundamental(luna, longitude, t) ≈ ERFA.faf03(t) atol=atol
        @test fundamental(luna, elongation, t) ≈ ERFA.fad03(t) atol=atol
        @test fundamental(luna, node, t) ≈ ERFA.faom03(t) atol=atol
        @test fundamental(mercury, t) ≈ ERFA.fame03(t) atol=atol
        @test fundamental(venus, t) ≈ ERFA.fave03(t) atol=atol
        @test fundamental(earth, t) ≈ ERFA.fae03(t) atol=atol
        @test fundamental(mars, t) ≈ ERFA.fama03(t) atol=atol
        @test fundamental(jupiter, t) ≈ ERFA.faju03(t) atol=atol
        @test fundamental(saturn, t) ≈ ERFA.fasa03(t) atol=atol
        @test fundamental(uranus, t) ≈ ERFA.faur03(t) atol=atol
        @test fundamental(neptune, t) ≈ ERFA.fane03(t) atol=atol
        @test fundamental(t) ≈ ERFA.fapa03(t) atol=atol
    end
    @testset "Obliquity" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))
        atol = 1e-12
        @test obliquity(iau1980, ep) ≈ ERFA.obl80(jd...) atol=atol
        @test obliquity(iau2006, ep) ≈ ERFA.obl06(jd...) atol=atol
    end
    @testset "Nutation" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))
        atol = 1e-12

        @testset "nut80" begin
            nut80_exp = ERFA.nut80(jd...)
            nut80_act = nutation(iau1980, ep)
            @test nut80_act[1] ≈ nut80_exp[1] atol=atol
            @test nut80_act[2] ≈ nut80_exp[2] atol=atol
        end
        @testset "nut00a" begin
            nut00a_exp = ERFA.nut00a(jd...)
            nut00a_act = nutation(iau2000a, ep)
            @test nut00a_act[1] ≈ nut00a_exp[1] atol=atol
            @test nut00a_act[2] ≈ nut00a_exp[2] atol=atol
        end
        @testset "nut00b" begin
            nut00b_exp = ERFA.nut00b(jd...)
            nut00b_act = nutation(iau2000b, ep)
            @test nut00b_act[1] ≈ nut00b_exp[1] atol=atol
            @test nut00b_act[2] ≈ nut00b_exp[2] atol=atol
        end
        @testset "nut06a" begin
            nut06a_exp = ERFA.nut06a(jd...)
            nut06a_act = nutation(iau2006a, ep)
            @test nut06a_act[1] ≈ nut06a_exp[1] atol=atol
            @test nut06a_act[2] ≈ nut06a_exp[2] atol=atol
        end
        @testset "numat" begin
            ϵ = obliquity(iau2006, ep)
            δψ, δϵ = nutation(iau2006a, ep)
            numat_exp = ERFA.numat(ϵ, δψ, δϵ)
            numat_act = nutation_matrix(ϵ, δψ, δϵ)
            @testset for i in eachindex(numat_act, numat_exp)
                @test numat_act[i] ≈ numat_exp[i] atol=atol
            end
        end
        @testset "nutm80" begin
            nutm80_exp = ERFA.nutm80(jd...)
            nutm80_act = nutation_matrix(iau1980, ep)
            @testset for i in eachindex(nutm80_act, nutm80_exp)
                @test nutm80_act[i] ≈ nutm80_exp[i] atol=atol
            end
        end
        @testset "num00a" begin
            num00a_exp = ERFA.num00a(jd...)
            num00a_act = nutation_matrix(iau2000a, ep)
            @testset for i in eachindex(num00a_act, num00a_exp)
                @test num00a_act[i] ≈ num00a_exp[i] atol=atol
            end
        end
        @testset "num00b" begin
            num00b_exp = ERFA.num00b(jd...)
            num00b_act = nutation_matrix(iau2000b, ep)
            @testset for i in eachindex(num00b_act, num00b_exp)
                @test num00b_act[i] ≈ num00b_exp[i] atol=atol
            end
        end
        @testset "num06a" begin
            num06a_exp = ERFA.num06a(jd...)
            num06a_act = nutation_matrix(iau2006a, ep)
            @testset for i in eachindex(num06a_act, num06a_exp)
                @test num06a_act[i] ≈ num06a_exp[i] atol=atol
            end
        end
    end
    @testset "Precession" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))
        atol = 1e-12

        @testset "bi00" begin
            bi00_exp = ERFA.bi00()
            bi00_act = bias(iau2000)
            @testset for i in eachindex(bi00_act, bi00_exp)
                @test bi00_act[i] ≈ bi00_exp[i] atol=atol
            end
        end
        @testset "pr00" begin
            pr00_exp = ERFA.pr00(jd...)
            pr00_act = precession(iau2000, ep)
            @test pr00_act[1] ≈ pr00_exp[1] atol=atol
            @test pr00_act[2] ≈ pr00_exp[2] atol=atol
        end
        @testset "bp00" begin
            bp00_exp = ERFA.bp00(jd...)
            bp00_act = bias_precession_matrix(iau2000, ep)
            @testset for i in eachindex(bp00_act)
                act = bp00_act[i]
                exp = bp00_exp[i]
                @testset for j in eachindex(act, exp)
                    @test act[j] ≈ exp[j] atol=atol
                end
            end
        end
        @testset "pfw06" begin
            pfw06_exp = ERFA.pfw06(jd...)
            pfw06_act = fukushima_williams(iau2006, ep)
            @testset for i in eachindex(pfw06_act)
                @test pfw06_act[i] ≈ pfw06_exp[i] atol=atol
            end
        end
        @testset "fw2m" begin
            pfw = fukushima_williams(iau2006, ep)
            fw2m_exp = ERFA.fw2m(pfw...)
            fw2m_act = fukushima_williams_matrix(pfw...)
            @testset for i in eachindex(fw2m_act)
                @test fw2m_act[i] ≈ fw2m_exp[i] atol=atol
            end
        end
    end
    @testset "Precession-Nutation" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        jd = value.(julian_twopart(TTEpoch(ep)))
        atol = 1e-12

        @testset "pn00" begin
            nut = nutation(iau2000a, ep)
            pn00_exp = ERFA.pn00(jd..., nut...)
            pn00_act = precession_nutation(iau2000, ep, nut...)
            @testset for i in eachindex(pn00_act, pn00_exp)
                if i == 1
                    @test pn00_act[i] ≈ pn00_exp[i] atol=atol
                else
                    act = pn00_act[i]
                    exp = pn00_exp[i]
                    @testset for j in eachindex(act, exp)
                        @test act[j] ≈ exp[j] atol=atol
                    end
                end
            end
        end
        @testset "pn00a" begin
            pn00a_exp = ERFA.pn00a(jd...)
            pn00a_act = precession_nutation(iau2000a, ep)
            @testset for i in eachindex(pn00a_act, pn00a_exp)
                if i in 1:3
                    @test pn00a_act[i] ≈ pn00a_exp[i] atol=atol
                else
                    act = pn00a_act[i]
                    exp = pn00a_exp[i]
                    @testset for j in eachindex(act, exp)
                        @test act[j] ≈ exp[j] atol=atol
                    end
                end
            end
        end
        @testset "pn00b" begin
            pn00b_exp = ERFA.pn00b(jd...)
            pn00b_act = precession_nutation(iau2000b, ep)
            @testset for i in eachindex(pn00b_act, pn00b_exp)
                if i in 1:3
                    @test pn00b_act[i] ≈ pn00b_exp[i] atol=atol
                else
                    act = pn00b_act[i]
                    exp = pn00b_exp[i]
                    @testset for j in eachindex(act, exp)
                        @test act[j] ≈ exp[j] atol=atol
                    end
                end
            end
        end
        @testset "pnm00a" begin
            pnm00a_exp = ERFA.pnm00a(jd...)
            pnm00a_act = precession_nutation_matrix(iau2000a, ep)
            @testset for i in eachindex(pnm00a_act, pnm00a_exp)
                @test pnm00a_act[i] ≈ pnm00a_exp[i] atol=atol
            end
        end
        @testset "pnm00b" begin
            pnm00b_exp = ERFA.pnm00b(jd...)
            pnm00b_act = precession_nutation_matrix(iau2000b, ep)
            @testset for i in eachindex(pnm00b_act, pnm00b_exp)
                @test pnm00b_act[i] ≈ pnm00b_exp[i] atol=atol
            end
        end
        @testset "pnm06a" begin
            pnm06a_exp = ERFA.pnm06a(jd...)
            pnm06a_act = precession_nutation_matrix(iau2006a, ep)
            @testset for i in eachindex(pnm06a_act, pnm06a_exp)
                @test pnm06a_act[i] ≈ pnm06a_exp[i] atol=atol
            end
        end
    end
    @testset "ICRS" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        tt = value.(julian_twopart(TTEpoch(ep)))
        ut = value.(julian_twopart(UT1Epoch(ep)))
        atol = 1e-12

        @testset "c2ixys" begin
            x, y, s = cip_coords_cio_locator(iau2000a, ep)
            c2ixys_act = celestial_to_intermediate(x, y, s)
            c2ixys_exp = ERFA.c2ixys(x, y, s)
            @testset for i in eachindex(c2ixys_act, c2ixys_exp)
                @test c2ixys_act[i] ≈ c2ixys_exp[i] atol=atol
            end
        end
        @testset "c2ixy" begin
            x, y, _ = cip_coords_cio_locator(iau2000a, ep)
            c2ixy_act = celestial_to_intermediate(iau2000, ep, x, y)
            c2ixy_exp = ERFA.c2ixy(tt..., x, y)
            @testset for i in eachindex(c2ixy_act, c2ixy_exp)
                @test c2ixy_act[i] ≈ c2ixy_exp[i] atol=atol
            end
        end
        @testset "c2ibpn" begin
            rnpb = precession_nutation_matrix(iau2000a, ep)
            c2ibpn_act = celestial_to_intermediate(iau2000, ep, rnpb)
            c2ibpn_exp = ERFA.c2ibpn(tt..., rnpb)
            @testset for i in eachindex(c2ibpn_act, c2ibpn_exp)
                @test c2ibpn_act[i] ≈ c2ibpn_exp[i] atol=atol
            end
        end
        @testset "c2i00a" begin
            c2i00a_act = celestial_to_intermediate(iau2000a, ep)
            c2i00a_exp = ERFA.c2i00a(tt...)
            @testset for i in eachindex(c2i00a_act, c2i00a_exp)
                @test c2i00a_act[i] ≈ c2i00a_exp[i] atol=atol
            end
        end
        @testset "c2i00b" begin
            c2i00b_act = celestial_to_intermediate(iau2000b, ep)
            c2i00b_exp = ERFA.c2i00b(tt...)
            @testset for i in eachindex(c2i00b_act, c2i00b_exp)
                @test c2i00b_act[i] ≈ c2i00b_exp[i] atol=atol
            end
        end
        @testset "c2tcio/c2teqx" begin
            xys = randn(3)
            a = polar_motion(iau2000, xys...)
            b = π
            c = polar_motion(iau2000, xys...)
            d_act = celestial_to_terrestrial_cio(a, b, c)
            e_act = celestial_to_terrestrial_equinox(a, b, c)
            d_exp = ERFA.c2tcio(a, b, c)
            e_exp = ERFA.c2teqx(a, b, c)
            @testset for i in eachindex(d_act, d_exp)
                @test d_act[i] ≈ d_exp[i] atol=atol
            end
            @testset for i in eachindex(e_act, e_exp)
                @test e_act[i] ≈ e_exp[i] atol=atol
            end
            @testset for i in eachindex(d_act, e_act)
                @test d_act[i] ≈ e_act[i] atol=atol
            end
        end
        @testset "c2t00a" begin
            xp = π/4
            yp = π/8
            c2t00a_act = celestial_to_terrestrial_cio(iau2000a, ep, xp, yp)
            c2t00a_exp = ERFA.c2t00a(tt..., ut..., xp, yp)
            @testset for i in eachindex(c2t00a_act, c2t00a_exp)
                @test c2t00a_act[i] ≈ c2t00a_exp[i] atol=atol
            end
        end
        @testset "c2t00b" begin
            xp = π/4
            yp = π/8
            c2t00b_act = celestial_to_terrestrial_cio(iau2000b, ep, xp, yp)
            c2t00b_exp = ERFA.c2t00b(tt..., ut..., xp, yp)
            @testset for i in eachindex(c2t00b_act, c2t00b_exp)
                @test c2t00b_act[i] ≈ c2t00b_exp[i] atol=atol
            end
        end
        @testset "c2tpe" begin
            xp = π/4
            yp = π/8
            δψ, δϵ = nutation(iau2000b, ep)
            c2tpe_act = celestial_to_terrestrial_equinox(iau2000, ep, δψ, δϵ, xp, yp)
            c2tpe_exp = ERFA.c2tpe(tt..., ut..., δψ, δϵ, xp, yp)
            @testset for i in eachindex(c2tpe_act, c2tpe_exp)
                @test c2tpe_act[i] ≈ c2tpe_exp[i] atol=atol
            end
        end
        @testset "c2txy" begin
            xp = float(π/4)
            yp = float(π/8)
            rnpb = precession_nutation_matrix(iau2000a, ep)
            x, y = cip_coords(rnpb)
            c2txy_act = celestial_to_terrestrial_cio(iau2000, ep, x, y, xp, yp)
            c2txy_exp = ERFA.c2txy(tt..., ut..., x, y, xp, yp)
            @testset for i in eachindex(c2txy_act, c2txy_exp)
                @test c2txy_act[i] ≈ c2txy_exp[i] atol=atol
            end
        end
        @testset "eors" begin
            rnpb = precession_nutation_matrix(iau2000a, ep)
            x, y = cip_coords(rnpb)
            s = cio_locator(iau2000, ep, x, y)
            @test equation_of_origins(rnpb, s) ≈ ERFA.eors(rnpb, s) atol=atol
        end
        @testset "pom00" begin
            xp = π/4
            yp = π/8
            s′ = tio_locator(iau2000, ep)
            pom00_act = polar_motion(iau2000, xp, yp, s′)
            pom00_exp = ERFA.pom00(xp, yp, s′)
            @testset for i in eachindex(pom00_act, pom00_exp)
                @test pom00_act[i] ≈ pom00_exp[i] atol=atol
            end
        end
        @testset "s00" begin
            rnpb = precession_nutation_matrix(iau2000a, ep)
            x, y = cip_coords(rnpb)
            @test cio_locator(iau2000, ep, x, y) ≈ ERFA.s00(tt..., x, y) atol=atol
        end
        @testset "s00a" begin
            @test cio_locator(iau2000a, ep) ≈ ERFA.s00a(tt...) atol=atol
        end
        @testset "s00b" begin
            @test cio_locator(iau2000b, ep) ≈ ERFA.s00b(tt...) atol=atol
        end
        @testset "s06" begin
            rnpb = precession_nutation_matrix(iau2000a, ep)
            x, y = cip_coords(rnpb)
            @test cio_locator(iau2006, ep, x, y) ≈ ERFA.s06(tt..., x, y) atol=atol
        end
        @testset "s06a" begin
            @test cio_locator(iau2006a, ep) ≈ ERFA.s06a(tt...) atol=atol
        end
        @testset "sp00" begin
            @test tio_locator(iau2000, ep) ≈ ERFA.sp00(tt...) atol=atol
        end
        @testset "xy06" begin
            x_act, y_act = cip_coords(iau2006, ep)
            x_exp, y_exp = ERFA.xy06(tt...)
            @test x_act ≈ x_exp atol=atol
            @test y_act ≈ y_exp atol=atol
        end
        @testset "xys00a" begin
            xys00a_act = cip_coords_cio_locator(iau2000a, ep)
            xys00a_exp = ERFA.xys00a(tt...)
            @testset for i in eachindex(xys00a_act, xys00a_exp)
                @test xys00a_act[i] ≈ xys00a_exp[i] atol=atol
            end
        end
        @testset "xys00b" begin
            xys00b_act = cip_coords_cio_locator(iau2000b, ep)
            xys00b_exp = ERFA.xys00b(tt...)
            @testset for i in eachindex(xys00b_act, xys00b_exp)
                @test xys00b_act[i] ≈ xys00b_exp[i] atol=atol
            end
        end
    end
    @testset "Rotation" begin
        ep = UTCEpoch(2020, 3, 16, 18, 15, 32.141)
        ut = value.(julian_twopart(UT1Epoch(ep)))
        tt = value.(julian_twopart(TTEpoch(ep)))
        tdb = value.(julian_twopart(TDBEpoch(ep)))
        atol = 1e-12

        @testset "era00" begin
            @test earth_rotation_angle(iau2000, ep) ≈ ERFA.era00(ut...) atol=atol
        end
        @testset "eqeq94" begin
            @test equinoxes(iau1994, ep) ≈ ERFA.eqeq94(tdb...) atol=atol
        end
        @testset "eect00" begin
            @test equinoxes(iau2000, ep) ≈ ERFA.eect00(tt...) atol=atol
        end
        @testset "ee00" begin
            ϵ = obliquity(iau1980, ep)
            δψ, _ = nutation(iau2000a, ep)
            @test equinoxes(iau2000, ep, ϵ, δψ) ≈ ERFA.ee00(tt..., ϵ, δψ) atol=atol
        end
        @testset "ee00a" begin
            @test equinoxes(iau2000a, ep) ≈ ERFA.ee00a(tt...) atol=atol
        end
        @testset "ee00b" begin
            @test equinoxes(iau2000b, ep) ≈ ERFA.ee00b(tt...) atol=atol
        end
        @testset "ee06a" begin
            @test equinoxes(iau2006a, ep) ≈ ERFA.ee06a(tt...) atol=atol
        end
        @testset "gmst82" begin
            @test mean_sidereal(iau1982, ep) ≈ ERFA.gmst82(ut...) atol=atol
        end
        @testset "gmst00" begin
            @test mean_sidereal(iau2000, ep) ≈ ERFA.gmst00(ut..., tt...) atol=atol
        end
        @testset "gmst06" begin
            @test mean_sidereal(iau2006, ep) ≈ ERFA.gmst06(ut..., tt...) atol=atol
        end
        @testset "gst94" begin
            @test apparent_sidereal(iau1994, ep) ≈ ERFA.gst94(ut...) atol=atol
        end
        @testset "gst00a" begin
            @test apparent_sidereal(iau2000a, ep) ≈ ERFA.gst00a(ut..., tt...) atol=atol
        end
        @testset "gst00b" begin
            @test apparent_sidereal(iau2000b, ep) ≈ ERFA.gst00b(ut...) atol=atol
        end
        @testset "gst06" begin
            rnpb = precession_nutation_matrix(iau2006a, ep)
            @test apparent_sidereal(iau2006, ep, rnpb) ≈ ERFA.gst06(ut..., tt..., rnpb) atol=atol
        end
        @testset "gst06a" begin
            @test apparent_sidereal(iau2006a, ep) ≈ ERFA.gst06a(ut..., tt...) atol=atol
        end
    end
end

