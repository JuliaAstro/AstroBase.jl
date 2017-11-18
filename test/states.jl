using JPLEphemeris
using RemoteFiles

de430 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp"
download(de430)
load_ephemeris!(SPK, path(de430))

@testset "States" begin
    r = [
    6068279.27,
    -1692843.94,
    -2516619.18,
    ]/1000

    v = [
    -660.415582,
    5495.938726,
    -5303.093233,
    ]/1000

    rv = [r; v]

    ep = UTCEpoch("2016-05-30T12:00:00.000")
    tdb = TDBEpoch(ep)

    iss = State(ep, r, v)
    rot = Rotation(GCRF, IAUEarth, ep)
    earth = state(ep, Earth)
    mars = state(ep, MarsBarycenter)
    mars_rot = State(tdb, r, v, IAUMars, MarsBarycenter)
    rv_mars = rv .+ earth .- mars
    r_mars = rv_mars[1:3]
    v_mars = rv_mars[4:6]

    @test iss == State(ep, r, v, GCRF, Earth)
    @test iss ≈ State(ep, r + eps(), v + eps(), GCRF, Earth)
    @test State(iss, frame=IAUEarth) ≈ State(ep, rot(r, v)..., IAUEarth)
    @test iss ≈ State(State(iss, frame=IAUEarth), frame=GCRF)
    @test State(iss, timescale=TDB) == State(tdb, r, v)
    @test State(iss, body=MarsBarycenter) ==
        State(ep, r_mars, v_mars, GCRF, MarsBarycenter)
    @test State(iss, timescale=TDB, body=MarsBarycenter) ==
        State(tdb, r_mars, v_mars, GCRF, MarsBarycenter)
    @test State(iss, frame=IAUEarth, body=MarsBarycenter) ≈
        State(ep, rot(r_mars, v_mars)..., IAUEarth, MarsBarycenter)
    @test State(iss, frame=IAUEarth, timescale=TDB) ≈
        State(tdb, rot(r, v)..., IAUEarth)
    @test State(State(mars_rot, body=Earth), body=MarsBarycenter) ≈ mars_rot
    @test State(State(mars_rot, frame=ITRF, body=Earth), frame=IAUMars,
        body=MarsBarycenter) ≈ mars_rot
    @test State(State(mars_rot, frame=ITRF, timescale=UT1, body=Earth),
        frame=IAUMars, timescale=TDB, body=MarsBarycenter) ≈ mars_rot
end
