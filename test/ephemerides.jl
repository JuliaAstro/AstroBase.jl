using AstroTime: TDBEpoch
using AstroBase.Ephemerides
using AstroBase.Bodies:
    NAIFId, ssb, sun, mercury, venus, earth, mars,
    jupiter, saturn, uranus, neptune

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
    # @testset "AbstractEphemeris" begin
    #     ep = TDBEpoch(2000, 1, 1)
    #     res = [1.0, 2.0, 3.0]
    #     @test position!(zeros(3), eph, ep, ssb, sun) == res
    #     @test position(eph, ep, ssb, sun) == res
    #     @test velocity!(zeros(3), eph, ep, ssb, sun) == res
    #     @test velocity(eph, ep, ssb, sun) == res
    #     @test state!(zeros(3), zeros(3), eph, ep, ssb, sun) == (res, res)
    #     @test state(eph, ep, ssb, sun) == (res, res)
    #
    #     @test position!(zeros(3), eph, ep, luna, io) == 4res
    #     @test position(eph, ep, luna, io) == 4res
    #     @test velocity!(zeros(3), eph, ep, luna, io) == 4res
    #     @test velocity(eph, ep, luna, io) == 4res
    #     @test state!(zeros(3), zeros(3), eph, ep, luna, io) == (4res, 4res)
    #     @test state(eph, ep, luna, io) == (4res, 4res)
    #
    #     @test position!(zeros(3), eph, ep, sun) == res
    #     @test position(eph, ep, sun) == res
    #     @test velocity!(zeros(3), eph, ep, sun) == res
    #     @test velocity(eph, ep, sun) == res
    #     @test state!(zeros(3), zeros(3), eph, ep, sun) == (res, res)
    #     @test state(eph, ep, sun) == (res, res)
    #
    #     @test position!(zeros(3), eph, ep, ssb) == zeros(3)
    #     @test position(eph, ep, ssb) == zeros(3)
    #     @test velocity!(zeros(3), eph, ep, ssb) == zeros(3)
    #     @test velocity(eph, ep, ssb) == zeros(3)
    #     @test state!(zeros(3), zeros(3), eph, ep, ssb) == (zeros(3), zeros(3))
    #     @test state(eph, ep, ssb) == (zeros(3), zeros(3))
    #
    #     @test position!(zeros(3), eph2, ep, luna, io) == 4res
    #     @test position(eph2, ep, luna, io) == 4res
    #     @test velocity!(zeros(3), eph2, ep, luna, io) == 4res
    #     @test velocity(eph2, ep, luna, io) == 4res
    #     @test state!(zeros(3), zeros(3), eph2, ep, luna, io) == (4res, 4res)
    #     @test state(eph2, ep, luna, io) == (4res, 4res)
    # end
    @testset "VSOP87" begin
        exp = Dict(sun => [-0.0017833873905459886,
                           0.007598638358181509,
                           -3.0846458879164256e-5,
                           -8.359521036032762e-6,
                           6.427858695952951e-7,
                           2.1628303284104923e-7],
                   mercury => [0.353046448850587,
                               -0.12833392955408252,
                               -0.04368943587287634,
                               0.0045755188619331704,
                               0.02754249508741387,
                               0.0018302598027261554],
                   venus => [0.6801118410706943,
                             -0.24261016437814215,
                             -0.04281468257088451,
                             0.006849335412957083,
                             0.018901335747826604,
                             -0.00013618119167514115],
                   earth => [-0.7157725461999801,
                             -0.7048878115760294,
                             5.6506882502700175e-6,
                             0.011870163513181561,
                             -0.012246037694124148,
                             3.246628790192791e-7],
                   mars => [-0.4459001810200729,
                            1.5510162330910005,
                            0.04320663825411424,
                            -0.01292770136907642,
                            -0.002679463170709558,
                            0.00026105673516857707],
                   jupiter => [-1.253636647162714,
                               -5.153732358252407,
                               0.04941516411266447,
                               0.007241201326338906,
                               -0.001423582458457613,
                               -0.0001560752464691477],
                   saturn => [2.5994179169904745,
                              -9.703657519703782,
                              0.06524583235735472,
                              0.005080974121526256,
                              0.0014262737807058702,
                              -0.00022706994269089197],
                   uranus => [16.75306441793305,
                              10.64563555849128,
                              -0.17750314479080884,
                              -0.002138293932829352,
                              0.003136269656333431,
                              3.9350166396732345e-5],
                   neptune => [29.07457968017663,
                               -7.09810729384933,
                               -0.52387424683584,
                               0.000723652561041086,
                               0.0030682008673629272,
                               -7.985981619277685e-5])

        ep = TDBEpoch(2019, 5, 6)
        vsop87 = VSOP87()
        @testset for (body, s_exp) in exp
            r_exp = s_exp[1:3]
            v_exp = s_exp[4:6]
            rv_exp = (r_exp, v_exp)
            r_act = position!(zeros(3), vsop87, ep, ssb, body)
            v_act = velocity!(zeros(3), vsop87, ep, ssb, body)
            rv_act = state!(zeros(3), zeros(3), vsop87, ep, ssb, body)
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
end
