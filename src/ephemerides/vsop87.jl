using MuladdMacro
using Parameters: @with_kw
using StaticArrays: @SMatrix

import ..AU

import AstroTime: Epoch, TDBEpoch, j2000, centuries, value, SECONDS_PER_CENTURY
import ..Bodies:
    SolarSystemBarycenter,
    Sun,
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune

export VSOP87

include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_sun.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_mercury.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_venus.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_earth.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_mars.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_jupiter.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_saturn.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_uranus.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_neptune.jl"))

const VSOP87_FK5 = @SMatrix [1.0 0.000000440360 -0.000000190919;
							 -0.000000479966 0.917482137087 -0.397776982902;
							 0.0 0.397776982902 0.917482137087]

@with_kw struct VSOP87 <: AbstractEphemeris
    order::Int = 5; @assert order >= 0; @assert order <= 5
    terms::Int = typemax(Int); @assert terms > 0
end

@inbounds function _vsop(coeffs, n_terms, max_order, max_terms, ep)
    t0 = value(centuries(j2000(TDBEpoch(ep)))) / 10.0

    # Powers of `t0`
    t = Array{Float64}(undef, 7)
    t[1] = 0.0 # t[-1]
    t[2] = 1.0 # t[0]
    t[3] = t0 # t[1]
    for i = 4:max_order + 2
        t[i] = t[i-1] * t0
    end

    r = zeros(3)
    v = zeros(3)

    for k = 1:3
        m = min(length(coeffs[k]) - 1, max_order)
        for order = 0:m
            i = order + 1
            it = order + 2
            n = min(n_terms[k][i], max_terms)
            for j = 1:n
                a = coeffs[k][i][1][j]
                b = coeffs[k][i][2][j]
                c = coeffs[k][i][3][j]

                u = b + c * t0
                su, cu = sincos(u)
                r[k] += a * cu * t[it]
                v[k] += order * a * cu * t[it-1] - a * c * su * t[it]
            end
        end
    end

    # AU to km
    r .*= AU
    # AU/millenium to km/s
    v .*= AU / 10SECONDS_PER_CENTURY

    VSOP87_FK5 * r, VSOP87_FK5 * v
end

for body in ("Sun", "Mercury", "Venus", "Earth", "Mars",
             "Jupiter", "Saturn", "Uranus", "Neptune")
    b = Symbol(body)
    c = Symbol("VSOP_", uppercase(body))
    n = Symbol("VSOP_", uppercase(body), "_NUM")
    @eval begin
        function position!(arr, eph::VSOP87, ep::Epoch, ::SolarSystemBarycenter, ::$b)
            _vsop($c, $n, eph.order, eph.terms, ep)[1]
        end

        function velocity!(arr, eph::VSOP87, ep::Epoch, ::SolarSystemBarycenter, ::$b)
            _vsop($c, $n, eph.order, eph.terms, ep)[2]
        end

        function state!(r0, v0, eph::VSOP87, ep::Epoch, ::SolarSystemBarycenter, ::$b)
            r1, v1 = _vsop($c, $n, eph.order, eph.terms, ep)
            r0 .+= r1
            v0 .+= v1
            r0, v0
        end
    end
end

function state!(pos, vel, vsop::VSOP87, ep::Epoch, from::SolarSystemBarycenter, to::CelestialBody)
    throw(ArgumentError("Body '$to' is not supported by VSOP87."))
end

function state!(pos, vel, vsop::VSOP87, ep::Epoch, from::CelestialBody, to::CelestialBody)
    rv1 = state!(zeros(3), zeros(3), vsop, ep, ssb, from)
    rv2 = state!(zeros(3), zeros(3), vsop, ep, ssb, to)
    pos .+= (rv2[1] .- rv1[1])
    vel .+= (rv2[2] .- rv1[2])
    pos, vel
end

function position!(pos, vsop::VSOP87, ep::Epoch, from::CelestialBody, to::CelestialBody)
    state!(pos, zeros(3), vsop, ep, from, to)
    pos
end

function velocity!(vel, vsop::VSOP87, ep::Epoch, from::CelestialBody, to::CelestialBody)
    state!(zeros(3), vel, vsop, ep, from, to)
    vel
end

