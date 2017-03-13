struct A <: Frame end
struct B <: Frame end
struct C <: Frame end
struct D <: Frame end

Rotation(::Type{A}, ::Type{B}, t::Float64) = Rotation{A,B}(diagm(1:6))
Rotation(::Type{B}, ::Type{C}, t::Float64) = Rotation{B,C}(diagm(1:6))
Rotation(::Type{C}, ::Type{D}, t::Float64) = Rotation{C,D}(diagm(1:6))
println(methods(Rotation))
# rotAD = Rotation(A, C, 0.0)
