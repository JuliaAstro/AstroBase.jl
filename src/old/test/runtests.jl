using AstroDynBase
using Test

AstroDynBase.update()

@testset "AstroDynBase" begin
    include("elements.jl")
    # include("rotations.jl")
    include("rotation_matrices.jl")
    include("bodies.jl")
    include("kepler.jl")
    #= include("ephemerides.jl") =#
end
