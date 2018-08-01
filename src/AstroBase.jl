module AstroBase

using Reexport

include("util.jl")
include(joinpath("bodies", "Bodies.jl"))
include("EarthAttitude.jl")

@reexport using .Bodies
@reexport using .EarthAttitude

export nutation_matrix80

"""
    nutation_matrix80(jd1, jd2)

Returns nutation matrix for a given 2 part Julian date (TDB).

# Example

```jldoctest
julia> nutation_matrix80(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotXZX{Float64}(0.409054, 3.75653e-5, -0.409017):
 1.0         -3.44666e-5  -1.494e-5
 3.44661e-5   1.0         -3.66564e-5
 1.49413e-5   3.66559e-5   1.0
 ```
"""
function nutation_matrix80(jd1, jd2)
    dpsi, deps = nutation(jd1, jd2)
    epsa = mean_obliquity_of_ecliptic(jd1, jd2)
    numat(epsa, dpsi, deps)
end
end # module
