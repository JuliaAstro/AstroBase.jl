module Ephemerides

using MuladdMacro

import AstroTime: Epoch, julian, centuries, value
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

function position!(arr, ::VSOP87, ep::Epoch, ::SolarSystemBarycenter, ::Sun; max_terms::Int=0)
    t = value(centuries(julian(ep))) / 10

    x = zeros(6)
    y = zeros(6)
    z = zeros(6)

    for i = 1:6
        for j in eachindex(sun_x[i][1])
            max_terms != 0 && j > max_terms && break

            @muladd x[i] += sun_x[i][1][j] * cos(sun_x[i][2][j] + sun_x[i][3][j] * t)
        end

        for j in eachindex(sun_y[i])
            max_terms != 0 && j > max_terms && break

            @muladd y[i] += sun_y[i][1][j] * cos(sun_y[i][2][j] + sun_y[i][3][j] * t)
        end

        for j in eachindex(sun_z[i])
            max_terms != 0 && j > max_terms && break

            @muladd z[i] += sun_x[i][1][j] * cos(sun_z[i][2][j] + sun_z[i][3][j] * t)
        end
    end

    arr[1] += @evalpoly t x[1] x[2] x[3] x[4] x[5] x[6]
    arr[2] += @evalpoly t y[1] y[2] y[3] y[4] y[5] y[6]
    arr[3] += @evalpoly t z[1] z[2] z[3] z[4] z[5] z[6]

    arr
end

end
