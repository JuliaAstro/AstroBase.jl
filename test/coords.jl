using AstroBase
using Test

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

    eph = VSOP87()

    ep = UTCEpoch("2016-05-30T12:00:00.000")
    tdb = TDBEpoch(ep)

    iss = State(ep, r, v)
    rot = Rotation(icrf, iau_earth, ep)
    pos_mars, vel_mars = state(eph, ep, earth, mars)
    r_mars = r .- pos_mars
    v_mars = v .- vel_mars
    mars_rot = State(tdb, r, v, frame=iau_mars, body=mars)

    @testset "Operators" begin
        @test iss == State(ep, r, v, frame=icrf, body=earth)
        @test iss ≈ State(ep, r .+ eps(), v .+ eps(), frame=icrf, body=earth)
    end
    @testset "Change frame" begin
        @test State(iss, eph, frame=iau_earth) ≈ State(ep, rot(r, v)..., frame=iau_earth)
        @test iss ≈ State(State(iss, eph, frame=iau_earth), eph, frame=icrf)
    end
    @testset "Change time scale" begin
        @test State(iss, eph, scale=TDB) == State(tdb, r, v)
    end
    @testset "Change central body" begin
        @test State(iss, eph, body=mars) ==
            State(ep, r_mars, v_mars, frame=icrf, body=mars)
        @test State(State(mars_rot, eph, body=earth), eph, body=mars) ≈ mars_rot
    end
    @testset "Change scale and body" begin
        @test State(iss, eph, scale=TDB, body=mars) ==
            State(tdb, r_mars, v_mars, frame=icrf, body=mars)
    end
    @testset "Change frame and body" begin
        @test State(iss, eph, frame=iau_earth, body=mars) ≈
            State(ep, rot(r_mars, v_mars)..., frame=iau_earth, body=mars)
        @test State(State(mars_rot, eph, frame=itrf, body=earth), eph, frame=iau_mars,
            body=mars) ≈ mars_rot
    end
    @testset "Change frame and scale" begin
        @test State(iss, eph, frame=iau_earth, scale=TDB) ≈
            State(tdb, rot(r, v)..., frame=iau_earth)
    end
    @testset "Change scale, frame, and body" begin
        @test State(State(mars_rot, eph, frame=itrf, scale=UT1, body=earth), eph,
            frame=iau_mars, scale=TDB, body=mars) ≈ mars_rot
    end
    @testset "Time Series" begin
        t = 1seconds:4seconds
        u = collect(1.0:4.0)
        ep = UTCEpoch(2000, 1, 1)
        ts = TimeSeries(ep, t, u)
        @test collect(ts) == [1.0, 2.0, 3.0, 4.0]
        @test ts[1] == 1.0
        @test ts(1.5seconds) ≈ 1.5
        @test ts(ep + 1.5seconds) ≈ 1.5
    end
    @testset "Trajectories" begin
        v_pos = [
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
        ]
        v_posvel = [
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
        ]
        v_posvelmass = [
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
        ]
        m_pos = [
            1.0 1.0 1.0;
            2.0 2.0 2.0;
            3.0 3.0 3.0;
        ]
        m_posvel = [
            1.0 1.0 1.0 1.0 1.0 1.0;
            2.0 2.0 2.0 2.0 2.0 2.0;
            3.0 3.0 3.0 3.0 3.0 3.0;
        ]
        m_posvelmass = [
            1.0 2.0 3.0;
            1.0 2.0 3.0;
            1.0 2.0 3.0;
            1.0 2.0 3.0;
            1.0 2.0 3.0;
            1.0 2.0 3.0;
        ]
    end
end

