using Test
import AstroBase

using AstroTime: Epoch, TDBEpoch, SECONDS_PER_DAY, j2000, seconds, value
using AstroBase.Ephemerides
using AstroBase.Bodies
using AstroBase: AU
using SPICE: spkezr, furnsh, kclear
using RemoteFiles: @RemoteFile, @RemoteFileSet, download, path

@RemoteFileSet KERNELS "SPK Kernels" begin
    de200 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de200.bsp"
    de245 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de245_1990_2010.bsp"
end

download(KERNELS)

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
    @testset "Meeus" begin
        furnsh(path(KERNELS, :de245))
        bodies = (
            mercury,
            venus,
            earth,
            mars,
            jupiter,
            saturn,
            uranus,
            neptune,
            pluto,
        )
        ep = TDBEpoch(2009, 5, 24)
        et = getet(ep)
        @testset for body in bodies
            name = string(body)
            if body in (jupiter, saturn, uranus, neptune, pluto)
                name *= " Barycenter"
            end
            s_exp = spkezr(name, et, "J2000", "NONE", "SUN")[1] .* 1e3
            r_exp = s_exp[1:3]
            v_exp = s_exp[4:6]
            rv_exp = (r_exp, v_exp)
            r_act = position!(zeros(3), meeus, ep, sun, body)
            v_act = velocity!(zeros(3), meeus, ep, sun, body)
            rv_act = state!(zeros(3), zeros(3), meeus, ep, sun, body)
            # TODO: Check tolerances
            @testset for i in 1:3
                @test r_act[i] ≈ r_exp[i] rtol=1e0
            end
            @testset for i in 1:3
                @test v_act[i] ≈ v_exp[i] rtol=1e-1
            end
            @testset for i in 1:3
                @test rv_act[1][i] ≈ rv_exp[1][i] rtol=1e0
                @test rv_act[2][i] ≈ rv_exp[2][i] rtol=1e-1
            end
        end
        kclear()
    end
    @testset "VSOP87" begin
        furnsh(path(KERNELS, :de200))
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
                @test r_act[i] ≈ r_exp[i] rtol=1e-2
            end
            @testset for i in 1:3
                @test v_act[i] ≈ v_exp[i] rtol=1e-3
            end
            @testset for i in 1:3
                @test rv_act[1][i] ≈ rv_exp[1][i] rtol=1e-2
                @test rv_act[2][i] ≈ rv_exp[2][i] rtol=1e-3
            end
        end
        kclear()
    end
end
