#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
using LinearAlgebra: I
using ReferenceFrameRotations: ddcm
using StaticArrays: SMatrix, SDiagonal

using ..Interfaces: AbstractFrame
using ..Time: Epoch

import Base: ∘

export
    Rotation,
    ∘,
    apply,
    apply!,
    compose,
    origin,
    target

struct Rotation{F1<:AbstractFrame,F2<:AbstractFrame,T,DT}
    origin::F1
    target::F2
    m::SMatrix{3,3,T}
    dm::SMatrix{3,3,DT}
end

function Rotation(origin, target, m::AbstractMatrix, dm::AbstractMatrix=zeros(3, 3))
    Rotation(origin, target, m, dm)
end

function Rotation(origin, target, m::AbstractMatrix, angular::AbstractVector)
    dm = ddcm(m, angular)
    return Rotation(origin, target, m, dm)
end

function Rotation(origin::F, target::F, ::Epoch) where F <: AbstractFrame
    Rotation(origin, target, Array{Float64}(I, 3, 3), zeros(3, 3))
end

origin(rot::Rotation) = rot.origin
target(rot::Rotation) = rot.target

function Base.inv(rot::Rotation)
    org = origin(rot)
    trg = target(rot)
    return Rotation(trg, org, rot.m', rot.dm')
end

function apply!(rot::Rotation, pos, vel)
    pos = rot.m * pos
    vel = rot.dm * pos + rot.m * vel
    return pos, vel
end

apply(rot::Rotation, pos, vel) = apply!(rot, copy(pos), copy(vel))

function apply!(rot::Rotation, vec)
    n = length(vec)
    n < 3 && throw(ArgumentError("`vec` must have at least 3 elements."))
    pos = vec[1:3]
    vec[1:3] = rot.m * pos
    n < 6 && return vec
    vec[4:6] = rot.dm * pos + rot.m * vec[4:6]
    return vec
end

apply(rot::Rotation, vec) = apply!(rot, copy(vec))

(rot::Rotation)(pos, vel) = apply(rot, pos, vel)
(rot::Rotation)(vec) = apply(rot, vec)

function compose(rot1::Rotation{F1, F}, rot2::Rotation{F, F2}) where {F, F1, F2}
    Rotation{F1, F2}(rot2.m * rot1.m, rot2.dm * rot1.m + rot2.m * rot1.dm)
end

∘(rot1::Rotation, rot2::Rotation) = compose(rot1, rot2)

function Rotation(frame1, frame2, ep::Epoch)
    frames = path_frames(frame1, frame2)
    rot = Rotation(frames[1], frames[2], ep)
    length(frames) == 2 && return rot

    for (f1, f2) in zip(frames[2:end-1], frames[3:end])
        rot = rot ∘ Rotation(from_sym(f1), from_sym(f2), ep)
    end
    rot
end

