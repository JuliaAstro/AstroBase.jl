#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

using LinearAlgebra: I
using ReferenceFrameRotations: ddcm
using StaticArrays: SMatrix, SVector

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
    function Rotation(
        origin::F1,
        target::F2,
        m::AbstractMatrix{T},
        dm::AbstractMatrix{DT}=zeros(T, 3, 3),
    ) where {F1,F2,T,DT}
        new{F1,F2,T,DT}(origin, target, m, dm)
    end
end

function Rotation(origin, target, m::AbstractMatrix, angular)
    dm = ddcm(m, SVector(angular))
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
    vel = rot.dm * pos + rot.m * vel
    pos = rot.m * pos
    return pos, vel
end

apply(rot::Rotation, pos, vel) = apply!(rot, copy(pos), copy(vel))
apply(rot::Rotation, posvel::Tuple) = apply(rot, posvel[1], posvel[2])

function apply!(rot::Rotation, vec::AbstractVector)
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
(rot::Rotation)(posvel::Tuple) = apply(rot, posvel[1], posvel[2])
(rot::Rotation)(vec::AbstractVector) = apply(rot, vec)

function compose(rot1::Rotation{F1, F}, rot2::Rotation{F, F2}) where {F, F1, F2}
    Rotation(F1(), F2(), rot2.m * rot1.m, rot2.dm * rot1.m + rot2.m * rot1.dm)
end

∘(rot1::Rotation, rot2::Rotation) = compose(rot1, rot2)

function Rotation(from::F1, to::F2, ep::Epoch) where {F1<:AbstractFrame,F2<:AbstractFrame}
    frames = path_frames(from, to)
    return _Rotation(ep, frames...)
end

@generated function _Rotation(ep::Epoch, path::AbstractFrame...)
    f1 = path[1]
    f2 = path[2]
    expr = :(Rotation($f1(), $f2(), ep))
    for i in 2:length(path) - 1
        f1 = path[i]
        f2 = path[i+1]
        expr = :(compose($expr, Rotation($f1(), $f2(), ep)))
    end
    return quote
        Base.@_inline_meta
        $expr
    end
end

