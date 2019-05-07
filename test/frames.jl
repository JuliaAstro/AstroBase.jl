using AstroBase
using SPICE: sxform, tisbod

@testset "Frames" begin
    @testset "Rotations" begin
        ep = TDBEpoch(2019, 5, 6)
        rot = Rotation(icrf, icrf, ep)
        @test origin(rot) == icrf
        @test target(rot) == icrf

        rot′ = inv(rot)
        @test rot == rot′

        r = randn(3)
        v = randn(3)

        rv = (r, v)
        r1, v1 = rot(r, v)
        rv1 = rot(rv)
        @test r1 == r
        @test v1 == v
        @test rv1 == rv

        comp = rot ∘ rot
        r1, v1 = comp(r, v)
        rv1 = comp(rv)
        @test r1 == r
        @test v1 == v
        @test rv1 == rv

        exp_path = [:IAUEarth, :ICRF, :CIRF, :TIRF, :ITRF]
        @test AstroBase.Frames.path_frames(iau_earth, itrf) == exp_path
    end
    @testset "Composition" begin
        ep = TDBEpoch(2019, 5, 6)
        et = value(seconds(j2000(ep)))
        rot = Rotation(iau_earth, iau_io, ep)
        m_act = rot.m
        m′_act = rot.m′
        s_exp = sxform("IAU_EARTH", "IAU_IO", et)
        m_exp = s_exp[1:3, 1:3]
        m′_exp = s_exp[4:6, 1:3]
        @testset for i in eachindex(m_act, m_exp)
            @test m_act[i] ≈ m_exp[i] atol=1e-8
        end
        @testset for i in eachindex(m′_act, m′_exp)
            @test m′_act[i] ≈ m′_exp[i] atol=1e-6
        end
    end
    @testset "IAU" begin
        ep = TDBEpoch(2019, 5, 6)
        et = value(seconds(j2000(ep)))
        @testset for body_name in AstroBase.Bodies.ALL_NAMES
            body = @eval $(Symbol(body_name))()
            id = naifid(body)
            if bodfnd(id, "PM")
                frame = @eval $(Symbol("IAU", body_name))()
                rot = Rotation(icrf, frame, ep)
                m_act = rot.m
                m′_act = rot.m′
                s_exp = tisbod("J2000", id, et)
                m_exp = s_exp[1:3, 1:3]
                m′_exp = s_exp[4:6, 1:3]
                @testset for i in eachindex(m_act, m_exp)
                    @test m_act[i] ≈ m_exp[i] atol=1e-8
                end
                @testset for i in eachindex(m′_act, m′_exp)
                    @test m′_act[i] ≈ m′_exp[i] atol=1e-6
                end
            end
        end
    end
    @testset "IERS" begin
        # Reference values from Orekit (http://www.orekit.org)
        tdb = TDBEpoch(2013, 3, 18, 12, 0, 0.0)
        r = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03]
        v = [7.37289205e+00, 2.08223573e+00, 4.39999794e-01]

        icrf_cirf = Rotation(icrf, cirf, tdb)
        cirf_icrf = Rotation(cirf, icrf, tdb)
        ref_icrf_cirf = [0.9999991427733881 -1.144503732190422E-8 -0.0013093710276537447;
            4.871197292135239E-8 0.9999999995949649 2.8461675919057616E-5;
            0.001309371026797659 -2.8461715302997654E-5 0.9999991423683543]
        r_cirf = [852.1380066867667, -4137.052915716736, 5296.806764985966]
        v_cirf = [7.372309565810324, 2.08224861625441, 0.4495940103957848]
        @test icrf_cirf.m ≈ ref_icrf_cirf
        @test all(icrf_cirf(r, v) .≈ (r_cirf, v_cirf))
        @test all(cirf_icrf(icrf_cirf(r, v)...) .≈ (r, v))
        @test all((cirf_icrf ∘ icrf_cirf)(r, v) .≈ (r, v))

        cirf_tirf = Rotation(cirf, tirf, tdb)
        tirf_cirf = Rotation(tirf, cirf, tdb)
        ref_cirf_tirf = [0.9972630428201841 -0.0739352650974288 0.0;
            0.0739352650974288 0.9972630428201841 0.0;
            0.0 0.0 0.9999999999999998]
        r_tirf = [1155.6798454967793, -4062.7269296118056, 5296.806764985964]
        v_tirf = [6.901921544550409, 2.537349749040834, 0.44959401039578467]
        @test cirf_tirf.m ≈ ref_cirf_tirf
        @test all(cirf_tirf(r_cirf, v_cirf) .≈ (r_tirf, v_tirf))
        @test all(tirf_cirf(cirf_tirf(r_cirf, v_cirf)...) .≈ (r_cirf, v_cirf))
        @test all((tirf_cirf ∘ cirf_tirf)(r_cirf, v_cirf) .≈ (r_cirf, v_cirf))

        tirf_itrf = Rotation(tirf, itrf, tdb)
        itrf_tirf = Rotation(itrf, tirf, tdb)
        ref_tirf_itrf = [0.9999999999999787 -3.010092334926554E-11 2.0567742763960954E-7;
            3.046034574888695E-11 0.999999999998473 -1.7475053232215554E-6;
            -2.0567742758669397E-7 1.7475053232277836E-6 0.9999999999984519]
        r_itrf = [1155.6809350526369, -4062.7361857684173, 5296.799427643569]
        v_itrf = [6.9019216369452225, 2.537348963579269, 0.4495970248578132]
        @test tirf_itrf.m ≈ ref_tirf_itrf
        @test all(tirf_itrf(r_tirf, v_tirf) .≈ (r_itrf, v_itrf))
        @test all(itrf_tirf(tirf_itrf(r_tirf, v_tirf)...) .≈ (r_tirf, v_tirf))
        @test all((itrf_tirf ∘ tirf_itrf)(r_tirf, v_tirf) .≈ (r_tirf, v_tirf))

        icrf_itrf = icrf_cirf ∘ cirf_tirf ∘ tirf_itrf
        r_exp, v_exp = tirf_itrf(cirf_tirf(icrf_cirf(r, v)))
        r_act, v_act = icrf_itrf(r, v)
        @testset for i in eachindex(r_exp, r_act)
            @test r_exp[i] ≈ r_act[i]
        end
        @testset for i in eachindex(v_exp, v_act)
            @test v_exp[i] ≈ v_act[i]
        end
    end
end
