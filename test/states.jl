using ERFA
using Compat

@testset "States" begin
    @testset "State type" begin
        ep = TTEpoch(2000, 1, 1)
        tdb = TDBEpoch(ep)
        @test State(ep, ones(6)) == State(ep, ones(6), GCRF, Earth)
        @test State(ep, ones(6)) ≈ State(ep, ones(6)+eps(), GCRF, Earth)
        rv = [ones(3)*7e3; ones(3)]
        s = State(TTEpoch(2000,1,1), ones(3)*7e3, ones(3))
        M = rotation_matrix(IAU{Earth}, GCRF, tdb)
        earth = state(Earth, tdb)
        mars = state(Mars, tdb)
        @test State(s, frame=IAU{Earth}) == State(ep, M*rv, IAU{Earth})
        @test s ≈ State(State(s, frame=IAU{Earth}), frame=GCRF)
        @test State(s, timescale=TDB) == State(tdb, rv)
        @test State(s, body=Mars) == State(ep, rv+earth-mars, GCRF, Mars)
        @test State(s, timescale=TDB, body=Mars) == State(tdb, rv+earth-mars, GCRF, Mars)
        @test State(s, frame=IAU{Earth}, body=Mars) ≈ State(ep, M*(rv+earth-mars), IAU{Earth}, Mars)
        @test State(s, frame=IAU{Earth}, timescale=TDB) == State(tdb, M*rv, IAU{Earth})
        mars_rot = State(tdb, rv, IAU{Mars}, Mars)
        @test State(State(mars_rot, body=Earth), body=Mars) ≈ mars_rot
        @test State(State(mars_rot, frame=ITRF, body=Earth), frame=IAU{Mars}, body=Mars) ≈ mars_rot
        @test State(State(mars_rot, frame=ITRF, timescale=UT1, body=Earth), frame=IAU{Mars}, timescale=TDB, body=Mars) ≈ mars_rot
    end
end
