#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

import ..Ephemerides: position, velocity, state
import ..Time: timescale

export AbstractCoord,
    AbstractState,
    AbstractTrajectory,
    centralbody,
    epoch,
    keplerian,
    position,
    state,
    refframe,
    timescale,
    velocity

abstract type AbstractCoord{S,F,B,T,ET,N} <: AbstractArray{T,N} end
abstract type AbstractState{S,F,B,T,ET} <: AbstractCoord{S,F,B,T,ET,1} end
abstract type AbstractTrajectory{S,F,B,T,ET} <: AbstractCoord{S,F,B,T,ET,2} end

function epoch end

timescale(::AbstractCoord{S}) where {S} = S()
refframe(::AbstractCoord{S,F}) where {S,F} = F()
centralbody(::AbstractCoord{S,F,B}) where {S,F,B} = B()

Base.eltype(::AbstractCoord{S,F,B,T}) where {S,F,B,T} = T
epochtype(::AbstractCoord{S,F,B,T,ET}) where {S,F,B,T,ET} = ET

Base.size(::AbstractState) = (6,)
Base.length(::AbstractState) = 6
Base.IndexStyle(::Type{<:AbstractState}) = IndexLinear()

