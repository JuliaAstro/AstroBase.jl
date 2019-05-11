module Ephemerides

using MuladdMacro

import AstroTime: Epoch, j2000, centuries, value, DAYS_PER_CENTURY
import ..Interfaces: AbstractEphemeris, position, position!
import ..Bodies: SolarSystemBarycenter, Sun

export VSOP87, vsop87, position!

include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_sun.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_mercury.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_venus.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_earth.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_mars.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_jupiter.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_saturn.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_uranus.jl"))
include(joinpath(@__DIR__, "..", "..", "gen", "vsop87_neptune.jl"))

struct VSOP87 <: AbstractEphemeris end
const vsop87 = VSOP87()

@inbounds function _vsop(coeffs, t0, n_terms, max_terms, max_order)
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
        for order = 0:max_order
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

    v ./= 10DAYS_PER_CENTURY

    r, v
end

function position!(arr, ::VSOP87, ep::Epoch, ::SolarSystemBarycenter, ::Sun)
    t0 = value(centuries(TDEpoch(ep))) / 10.0
    _vsop(VSOP_SUN, t0, VSOP_SUN_NUM, typemax(Int), 5)[1]
end

end
