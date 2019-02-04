using AstroTime: TDBEpoch
using AstroBase.Interfaces
using AstroBase.Bodies: ssb, sun, earth, NAIFId

struct DummyEphemeris <: AbstractEphemeris end
const eph = DummyEphemeris()

function AstroBase.Interfaces.position!(pos, ::DummyEphemeris, ::Epoch, ::NAIFId, ::NAIFId)
    pos .+= [1.0, 2.0, 3.0]
end

function AstroBase.Interfaces.velocity!(vel, ::DummyEphemeris, ::Epoch, ::NAIFId, ::NAIFId)
    vel .+= [1.0, 2.0, 3.0]
end

@testset "Interfaces" begin
    @testset "AbstractEphemeris" begin
        ep = TDBEpoch(2000, 1, 1)
        res = [1.0, 2.0, 3.0]
        @test position!(zeros(3), eph, ep, ssb, sun) == res
        @test position(eph, ep, ssb, sun) == res
        @test velocity!(zeros(3), eph, ep, ssb, sun) == res
        @test velocity(eph, ep, ssb, sun) == res
        @test position_velocity!(zeros(3), zeros(3), eph, ep, ssb, sun) == (res, res)
        @test position_velocity(eph, ep, ssb, sun) == (res, res)

        @test position!(zeros(3), eph, ep, sun) == res
        @test position(eph, ep, sun) == res
        @test velocity!(zeros(3), eph, ep, sun) == res
        @test velocity(eph, ep, sun) == res
        @test position_velocity!(zeros(3), zeros(3), eph, ep, sun) == (res, res)
        @test position_velocity(eph, ep, sun) == (res, res)

        @test position!(zeros(3), eph, ep, ssb) == zeros(3)
        @test position(eph, ep, ssb) == zeros(3)
        @test velocity!(zeros(3), eph, ep, ssb) == zeros(3)
        @test velocity(eph, ep, ssb) == zeros(3)
        @test position_velocity!(zeros(3), zeros(3), eph, ep, ssb) == (zeros(3), zeros(3))
        @test position_velocity(eph, ep, ssb) == (zeros(3), zeros(3))
    end
end
