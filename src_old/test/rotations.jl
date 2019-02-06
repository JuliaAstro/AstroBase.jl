import AstroDynBase.Rotation

struct A <: Frame end
struct B <: Frame end
struct C <: Frame end
struct D <: Frame end

Rotation(::Type{A}, ::Type{B}, ep) = Rotation{A,B}(diagm(1:6))
Rotation(::Type{B}, ::Type{C}, ep) = Rotation{B,C}(diagm(1:6))
Rotation(::Type{C}, ::Type{D}, ep) = Rotation{C,D}(diagm(1:6))

rotAD = Rotation(A, D, 0.0)

@testset "Rotations" begin
    @test rotAD.matrix == diagm((1.0:6.0).^3)
end
