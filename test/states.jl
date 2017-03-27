@testset "States" begin
    r = [
    6068279.27,
    -1692843.94,
    -2516619.18,
    ]/1000*km

    v = [
    -660.415582,
    5495.938726,
    -5303.093233,
    ]/1000*kps

    ep = UTCEpoch("2016-05-30T12:00:00.000")
    tdb = TDBEpoch(ep)

    iss = State(ep, r, v)

    @test iss == State(ep, r, v, GCRF, Earth)
    @test iss ≈ State(ep, r + eps()*km, v + eps()*kps, GCRF, Earth)
    rot = Rotation(GCRF, IAUEarth, ep)
    #= earth = state(Earth, tdb) =#
    #= mars = state(Mars, tdb) =#
    @test State(iss, frame=IAUEarth) == State(ep, rot(r, v)..., IAUEarth)
    #= @test s ≈ State(State(s, frame=IAU{Earth}), frame=GCRF) =#
    @test State(iss, timescale=TDB) == State(tdb, r, v)
    #= @test State(s, body=Mars) == State(ep, rv+earth-mars, GCRF, Mars) =#
    #= @test State(s, timescale=TDB, body=Mars) == State(tdb, rv+earth-mars, GCRF, Mars) =#
    #= @test State(s, frame=IAU{Earth}, body=Mars) ≈ State(ep, M*(rv+earth-mars), IAU{Earth}, Mars) =#
    @test State(iss, frame=IAUEarth, timescale=TDB) == State(tdb, rot(r, v)..., IAUEarth)
    #= mars_rot = State(tdb, rv, IAU{Mars}, Mars) =#
    #= @test State(State(mars_rot, body=Earth), body=Mars) ≈ mars_rot =#
    #= @test State(State(mars_rot, frame=ITRF, body=Earth), frame=IAU{Mars}, body=Mars) ≈ mars_rot =#
    #= @test State(State(mars_rot, frame=ITRF, timescale=UT1, body=Earth), frame=IAU{Mars}, timescale=TDB, body=Mars) ≈ mars_rot =#
end
