using LinearAlgebra: I
using ReferenceFrameRotations: compose_rotation
using StaticArrays: SMatrix, SDiagonal

export Rotation

struct Rotation{T, TV}
    m::SMatrix{3, 3, T}
    δm::SMatrix{3, 3, TV}
end

Rotation(::F, ::F, ::Epoch) where {F} = Rotation{Float64, Float64}(Array{Float64}(I, 3, 3), Array{Float64}(I, 3, 3))
#=  =#
#= function Rotation(::Frame1, ::Frame2) where {Frame1, Frame2} =#
#=     frames = path_frames(Frame1, Frame2) =#
#= end =#
