module Coords

using AstroTime: Epoch, TimeScale

import ..Interfaces:
    position,
    position_velocity,
    velocity

using ..Bodies:
    CelestialBody,
    earth,
    grav_param

using ..Frames:
    AbstractFrame,
    icrf

import ..TwoBody:
    cartesian,
    keplerian

using StaticArrays: SVector

export State, KeplerianState, icrf, epoch, frame, body

#################
# AbstractState #
#################

export AbstractState, epoch, frame

abstract type AbstractState{Scale, Frame, Body} end

function epoch end

frame(s::AbstractState{_S, Frame}) where {_S, Frame} = Frame
body(s::AbstractState{_S, _F, Body}) where {_S, _F, Body} = Body

struct State{Scale, Frame, Body, T, TP, TV} <: AbstractState{Scale, Frame, Body}
    epoch::Epoch{Scale, T}
    pos::SVector{3, TP}
    vel::SVector{3, TV}
    function State(epoch::Epoch{Scale, T}, pos, vel;
                   frame::Frame=icrf, body::Body=earth) where {Scale, Frame, Body, T}
        TP = eltype(pos)
        TV = eltype(vel)
        new{Scale::TimeScale, frame::AbstractFrame, body::CelestialBody, T, TP, TV}(epoch, pos, vel)
    end
end

epoch(s::State) = s.epoch
position(s::State) = s.pos
velocity(s::State) = s.vel
position_velocity(s::State) = (s.pos, s.vel)
keplerian(s::State) = keplerian(position(s), velocity(s), grav_param(body(s)))

struct KeplerianState{Scale, Frame, Body, T} <: AbstractState{Scale, Frame, Body}
    epoch::Epoch{Scale, T}
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    ν::T
    function KeplerianState(epoch::Epoch{Scale, T}, a, e, i, Ω, ω, ν;
                            frame::Frame=icrf, body::Body=earth) where {Scale, Frame, Body, T}
        new{Scale::TimeScale, frame::AbstractFrame, body::CelestialBody, T}(epoch, a, e, i, Ω, ω, ν)
    end
end

epoch(s::KeplerianState) = s.epoch
position(s::KeplerianState) = cartesian(s.a, s.i, s.e, s.Ω, s.ω, s.ν, grav_param(body(s)))[1]
velocity(s::KeplerianState) = cartesian(s.a, s.i, s.e, s.Ω, s.ω, s.ν, grav_param(body(s)))[2]
position_velocity(s::KeplerianState) = cartesian(s.a, s.i, s.e, s.Ω, s.ω, s.ν, grav_param(body(s)))
keplerian(s::KeplerianState) = (s.a, s.i, s.e, s.Ω, s.ω, s.ν)

KeplerianState(s::AbstractState) = KeplerianState(epoch(s), keplerian(s)...; frame=frame(s), body=body(s))
State(s::KeplerianState) = State(epoch(s), position_velocity(s)...; frame=frame(s), body=body(s))

struct Trajectory{Scale, Frame, Body}
end

end
