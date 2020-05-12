using AstroBase
using LinearAlgebra: norm
using OrdinaryDiffEq: ODEProblem, solve, Vern9
using Test

function rhs(du, u, p, t)
    du[1:3] = u[4:6]
    r3 = norm(u[1:3])^3
    du[4:6] = -p .* u[1:3] ./ r3
end

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

ele = keplerian(r, v, grav_param(earth))
t1 = period(ele[1], grav_param(earth))
s0_exp = State(ep, r, v)
k0_exp = KeplerianState(s0_exp)

@testset "Coords" begin
    @testset "States" begin
        eph = VSOP87()

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
    # @testset "Trajectories" begin
    #     @testset "Transpose - Array of Arrays" begin
    #         a = [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]
    #         b = [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    #         @test AstroBase.Coords.transpose_aoa(a) == b
    #     end
    #     @testset "Propagation" begin
    #         prob = ODEProblem(rhs, rv, (0.0, value(t1)), p=grav_param(earth))
    #         sol = solve(prob, Vern9())
    #         tra = Trajectory(ep, sol.t * seconds, sol.u)
    #         rv0 = tra(0.0seconds)
    #         ele0 = keplerian(tra, 0.0seconds)
    #         s0 = State(tra, 0.0seconds)
    #         k0 = KeplerianState(tra, 0.0seconds)
    #         @test rv0 ≈ rv
    #         @test collect(ele0) ≈ collect(ele)
    #         @test s0_exp ≈ s0
    #         @test k0_exp ≈ k0
    #         rv1 = tra(t1)
    #         ele1 = keplerian(tra, t1)
    #         s1 = State(tra, t1)
    #         k1 = KeplerianState(tra, t1)
    #         @test rv1 ≈ rv atol=0.1
    #         @test collect(ele1) ≈ collect(ele) atol=0.1
    #         @test s0_exp ≈ s1 atol=0.1
    #         @test k0_exp ≈ k1 atol=0.1
    #     end
    # end
end

