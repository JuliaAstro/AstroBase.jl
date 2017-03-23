@testset "IERS" begin
    # Reference values from Orekit (http://www.orekit.org)
    tdb = TDBEpoch(2013, 3, 18, 12, 0, 0.0)
    rv = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03, 7.37289205e+00,
        2.08223573e+00, 4.39999794e-01]

    gcrf_cirf = Rotation(GCRF, CIRF, tdb)
    cirf_gcrf = Rotation(CIRF, GCRF, tdb)
    ref_gcrf_cirf = [0.9999991427733881 -1.144503732190422E-8 -0.0013093710276537447;
        4.871197292135239E-8 0.9999999995949649 2.8461675919057616E-5;
        0.001309371026797659 -2.8461715302997654E-5 0.9999991423683543]
    rv_cirf = [852.1380066867667, -4137.052915716736, 5296.806764985966,
        7.372309565810324, 2.08224861625441, 0.4495940103957848]
    @test gcrf_cirf isa Rotation{GCRF,CIRF}
    @test gcrf_cirf.matrix[1:3,1:3] ≈ ref_gcrf_cirf
    @test gcrf_cirf(rv) ≈ rv_cirf
    @test cirf_gcrf isa Rotation{CIRF,GCRF}
    @test cirf_gcrf(gcrf_cirf(rv)) ≈ rv

    cirf_tirf = Rotation(CIRF, TIRF, tdb)
    tirf_cirf = Rotation(TIRF, CIRF, tdb)
    ref_cirf_tirf = [0.9972630428201841 -0.0739352650974288 0.0;
        0.0739352650974288 0.9972630428201841 0.0;
        0.0 0.0 0.9999999999999998]
    rv_tirf = [1155.6798454967793, -4062.7269296118056, 5296.806764985964,
        6.901921544550409, 2.537349749040834, 0.44959401039578467]
    @test cirf_tirf isa Rotation{CIRF,TIRF}
    @test cirf_tirf.matrix[1:3,1:3] ≈ ref_cirf_tirf
    @test cirf_tirf(rv_cirf) ≈ rv_tirf
    @test tirf_cirf isa Rotation{TIRF,CIRF}
    @test tirf_cirf(cirf_tirf(rv_cirf)) ≈ rv_cirf

    tirf_itrf = Rotation(TIRF, ITRF, tdb)
    itrf_tirf = Rotation(ITRF, TIRF, tdb)
    ref_tirf_itrf = [0.9999999999999787 -3.010092334926554E-11 2.0567742763960954E-7;
        3.046034574888695E-11 0.999999999998473 -1.7475053232215554E-6;
        -2.0567742758669397E-7 1.7475053232277836E-6 0.9999999999984519]
    rv_itrf = [1155.6809350526369, -4062.7361857684173, 5296.799427643569,
        6.9019216369452225, 2.537348963579269, 0.4495970248578132]
    @test tirf_itrf isa Rotation{TIRF,ITRF}
    @test tirf_itrf.matrix[1:3,1:3] ≈ ref_tirf_itrf
    @test tirf_itrf(rv_tirf) ≈ rv_itrf
    @test itrf_tirf isa Rotation{ITRF,TIRF}
    @test itrf_tirf(tirf_itrf(rv_tirf)) ≈ rv_tirf

    #= s_gcrf = State(tdb, rv) =#
    #= s_cirf = State(tdb, rv_cirf, CIRF) =#
    #= s_tirf = State(tdb, rv_tirf, TIRF) =#
    #= s_itrf = State(tdb, rv_itrf, ITRF) =#
    #= @test s_gcrf ≈ State(State(s_gcrf, frame=CIRF), frame=GCRF) =#
    #= @test s_cirf ≈ State(State(s_cirf, frame=TIRF), frame=CIRF) =#
    #= @test s_tirf ≈ State(State(s_tirf, frame=ITRF), frame=TIRF) =#
    #= @test s_gcrf ≈ State(State(s_gcrf, frame=ITRF), frame=GCRF) =#
end
