using Test
import AstroBase

using AstroTime: Epoch, TDBEpoch, SECONDS_PER_DAY, j2000, seconds, value, julian_twopart
using AstroBase.Ephemerides
using AstroBase.Bodies
using AstroBase: AU
using ERFA
using LinearAlgebra: norm
using SPICE: spkezr, furnsh, kclear
using RemoteFiles: @RemoteFile, @RemoteFileSet, download, path

@RemoteFile(
    DE200,
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de200.bsp",
    dir=joinpath(@__DIR__(), "data"),
)

download(DE200)
furnsh(path(DE200))

getet(epoch) = value(seconds(j2000(epoch)))

import AstroBase.Ephemerides: position!, velocity!, state!

struct DummyEphemeris <: AbstractEphemeris end
const eph = DummyEphemeris()

function position!(pos, ::DummyEphemeris, ::Epoch, ::NAIFId, ::NAIFId)
    pos .+= [1.0, 2.0, 3.0]
end

function velocity!(vel, ::DummyEphemeris, ::Epoch, ::NAIFId, ::NAIFId)
    vel .+= [1.0, 2.0, 3.0]
end

function state!(pos, vel, ::DummyEphemeris, ::Epoch, ::NAIFId, ::NAIFId)
    (pos .+= [1.0, 2.0, 3.0], vel .+= [1.0, 2.0, 3.0])
end

struct DummyEphemeris2 <: AbstractEphemeris end
const eph2 = DummyEphemeris2()

function position!(pos, ::DummyEphemeris2, ::Epoch, ::NAIFId, ::NAIFId)
    pos .+= [1.0, 2.0, 3.0]
end

function velocity!(vel, ::DummyEphemeris2, ::Epoch, ::NAIFId, ::NAIFId)
    vel .+= [1.0, 2.0, 3.0]
end

@testset "Ephemerides" begin
    @testset "AbstractEphemeris" begin
        ep = TDBEpoch(2000, 1, 1)
        res = [1.0, 2.0, 3.0]
        @test position!(zeros(3), eph, ep, ssb, sun) == res
        @test position(eph, ep, ssb, sun) == res
        @test velocity!(zeros(3), eph, ep, ssb, sun) == res
        @test velocity(eph, ep, ssb, sun) == res
        @test state!(zeros(3), zeros(3), eph, ep, ssb, sun) == (res, res)
        @test state(eph, ep, ssb, sun) == (res, res)

        @test position!(zeros(3), eph, ep, luna, io) == 4res
        @test position(eph, ep, luna, io) == 4res
        @test velocity!(zeros(3), eph, ep, luna, io) == 4res
        @test velocity(eph, ep, luna, io) == 4res
        @test state!(zeros(3), zeros(3), eph, ep, luna, io) == (4res, 4res)
        @test state(eph, ep, luna, io) == (4res, 4res)

        @test position!(zeros(3), eph, ep, sun) == res
        @test position(eph, ep, sun) == res
        @test velocity!(zeros(3), eph, ep, sun) == res
        @test velocity(eph, ep, sun) == res
        @test state!(zeros(3), zeros(3), eph, ep, sun) == (res, res)
        @test state(eph, ep, sun) == (res, res)

        @test position!(zeros(3), eph, ep, ssb) == zeros(3)
        @test position(eph, ep, ssb) == zeros(3)
        @test velocity!(zeros(3), eph, ep, ssb) == zeros(3)
        @test velocity(eph, ep, ssb) == zeros(3)
        @test state!(zeros(3), zeros(3), eph, ep, ssb) == (zeros(3), zeros(3))
        @test state(eph, ep, ssb) == (zeros(3), zeros(3))

        @test position!(zeros(3), eph2, ep, luna, io) == 4res
        @test position(eph2, ep, luna, io) == 4res
        @test velocity!(zeros(3), eph2, ep, luna, io) == 4res
        @test velocity(eph2, ep, luna, io) == 4res
        @test state!(zeros(3), zeros(3), eph2, ep, luna, io) == (4res, 4res)
        @test state(eph2, ep, luna, io) == (4res, 4res)
    end
    @testset "Simon/Bretagnon" begin
        bodies = (
            mercury,
            venus,
            earth_barycenter,
            mars,
            jupiter,
            saturn,
            uranus,
            neptune,
        )
        pos_error = (
            0.5e6,
            1.1e6,
            1.3e6,
            9e6,
            82e6,
            263e6,
            661e6,
            248e6,
        )
        vel_error = (
            0.7,
            0.9,
            1.0,
            2.5,
            8.2,
            24.6,
            27.4,
            21.4,
        )
        ep = TDBEpoch(2009, 5, 24)
        et = getet(ep)
        jd1, jd2 = value.(julian_twopart(ep))
        @testset for (id, body) in enumerate(bodies)
            @testset "ERFA" begin
                r_erfa, v_erfa = ERFA.plan94(jd1, jd2, id)
                r_erfa .*= AU
                v_erfa .*= AU
                v_erfa ./= SECONDS_PER_DAY
                rv_erfa = (r_erfa, v_erfa)
                r_act = position!(zeros(3), simon_bretagnon, ep, sun, body)
                v_act = velocity!(zeros(3), simon_bretagnon, ep, sun, body)
                rv_act = state!(zeros(3), zeros(3), simon_bretagnon, ep, sun, body)
                @testset for i in 1:3
                    @test r_act[i] ≈ r_erfa[i] atol=1.0
                end
                @testset for i in 1:3
                    @test v_act[i] ≈ v_erfa[i] atol=5.0
                end
                @testset for i in 1:3
                    @test rv_act[1][i] ≈ rv_erfa[1][i] atol=1.0
                    @test rv_act[2][i] ≈ rv_erfa[2][i] atol=5.0
                end
            end
            @testset "JPL" begin
                name = string(body)
                if body in (jupiter, saturn, uranus, neptune)
                    name *= " Barycenter"
                end
                s_exp = spkezr(name, et, "J2000", "NONE", "SUN")[1] .* 1e3
                r_exp = s_exp[1:3]
                v_exp = s_exp[4:6]
                rv_exp = (r_exp, v_exp)
                r_act = position!(zeros(3), simon_bretagnon, ep, sun, body)
                v_act = velocity!(zeros(3), simon_bretagnon, ep, sun, body)
                rv_act = state!(zeros(3), zeros(3), simon_bretagnon, ep, sun, body)
                @test norm(r_act) ≈ norm(r_exp) atol=pos_error[id]
                @test norm(v_act) ≈ norm(v_exp) atol=vel_error[id]
                @test norm(rv_act[1]) ≈ norm(rv_exp[1]) atol=pos_error[id]
                @test norm(rv_act[2]) ≈ norm(rv_exp[2]) atol=vel_error[id]
            end
        end
    end
    @testset "VSOP87" begin
        bodies = (
            sun,
            mercury,
            venus,
            earth,
            mars,
            jupiter,
            saturn,
            uranus,
            neptune,
        )
        ep = TDBEpoch(2019, 5, 6)
        et = getet(ep)
        vsop87 = VSOP87()
        @testset for body in bodies
            name = string(body)
            if body in (jupiter, saturn, uranus, neptune)
                name *= " Barycenter"
            end
            s_exp = spkezr(name, et, "J2000", "NONE", "SSB")[1] .* 1e3
            r_exp = s_exp[1:3]
            v_exp = s_exp[4:6]
            rv_exp = (r_exp, v_exp)
            r_act = position!(zeros(3), vsop87, ep, ssb, body)
            v_act = velocity!(zeros(3), vsop87, ep, ssb, body)
            rv_act = state!(zeros(3), zeros(3), vsop87, ep, ssb, body)
            # TODO: Check tolerances
            @testset for i in 1:3
                @test r_act[i] ≈ r_exp[i]
            end
            @testset for i in 1:3
                @test v_act[i] ≈ v_exp[i]
            end
            @testset for i in 1:3
                @test rv_act[1][i] ≈ rv_exp[1][i]
                @test rv_act[2][i] ≈ rv_exp[2][i]
            end
        end
    end

    kclear()
end
