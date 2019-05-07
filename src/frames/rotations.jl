import Base: ∘

using AstroTime: Epoch

using LinearAlgebra: I
using ReferenceFrameRotations: compose_rotation
using StaticArrays: SMatrix, SDiagonal

export
    Rotation,
    ∘,
    compose,
    origin,
    target

struct Rotation{F1, F2, T, T′}
    m::SMatrix{3, 3, T}
    m′::SMatrix{3, 3, T′}
    function Rotation{F1, F2}(m, m′=zeros(3, 3)) where {F1, F2}
        T = eltype(m)
        T′ = eltype(m′)
        new{F1, F2, T, T′}(m, m′)
    end
end

function Rotation(origin::F, target::F, ::Epoch) where {F}
    Rotation{origin::AbstractFrame,
             target::AbstractFrame}(Array{Float64}(I, 3, 3), zeros(3, 3))
end

origin(::Rotation{F}) where {F} = F
target(::Rotation{_F, F}) where {_F, F} = F

function Base.inv(rot::Rotation{F1, F2}) where {F1, F2}
    Rotation{F2, F1}(rot.m', rot.m′')
end

(rot::Rotation)(r, v) = rot.m * r, rot.m′ * r + rot.m * v
(rot::Rotation)(rv) = rot.m * rv[1], rot.m′ * rv[1] + rot.m * rv[2]

function compose(rot1::Rotation{F1, F}, rot2::Rotation{F, F2}) where {F, F1, F2}
    Rotation{F1, F2}(rot2.m * rot1.m, rot2.m′ * rot1.m + rot2.m * rot1.m′)
end

∘(rot1::Rotation, rot2::Rotation) = compose(rot1, rot2)

function Rotation(frame1, frame2, ep::Epoch)
    frames = path_frames(frame1, frame2)
    rot = Rotation(from_sym(frames[1]), from_sym(frames[2]), ep)
    length(frames) == 2 && return rot

    for (f1, f2) in zip(frames[2:end-1], frames[3:end])
        rot = rot ∘ Rotation(from_sym(f1), from_sym(f2), ep)
    end
    rot
end
