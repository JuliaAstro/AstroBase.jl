using LinearAlgebra: promote_leaf_eltypes
using StaticArrays: SVector

using ..Bodies: CelestialBody, earth, grav_param
using ..Frames: AbstractFrame, icrf
using ..Time: Epoch, TimeScale

import Base: ==
import ..Interfaces: AbstractState,
    centralbody,
    epoch,
    position,
    refframe,
    state,
    timescale,
    velocity
import ..TwoBody: period, cartesian, keplerian

export State, KeplerianState, pos_type, vel_type, period

function period(s::AbstractState)
    μ = grav_param(body(s))
    a, _ = keplerian(s)
    return period(a, μ)
end

struct State{S,F,B,T,ET} <: AbstractState{S,F,B,T,ET}
    epoch::Epoch{S,ET}
    frame::F
    body::B
    pos::SVector{3,T}
    vel::SVector{3,T}

    function State(
        epoch::Epoch{S1,ET}, pos, vel;
        scale::S2=timescale(epoch),
        frame::F=icrf,
        body::B=earth
    ) where {S1,S2,F,B,ET}
        TP = eltype(pos)
        TV = eltype(vel)
        T = Union{TP,TV}
        ep = Epoch{S2}(epoch)
        new{S2,F,B,T,ET}(ep, frame, body, pos, vel)
    end
end

function Base.getindex(s::State, i::Int)
    return ifelse(i <= 3, s.pos[i], s.vel[i])
end

function (==)(s1::State, s2::State)
    s1.epoch == s2.epoch && s1.pos == s2.pos && s1.vel == s2.vel
end

function Base.isapprox(
    s1::State,
    s2::State,
    atol::Real=0,
    rtol::Real=Base.rtoldefault(
        promote_leaf_eltypes(s1),
        promote_leaf_eltypes(s2),
        atol,
    ),
)
    isapprox(s1.epoch, s2.epoch, atol=atol, rtol=rtol) &&
    isapprox(s1.pos, s2.pos, atol=atol, rtol=rtol) &&
    isapprox(s1.vel, s2.vel, atol=atol, rtol=rtol)
end

epoch(s::State) = s.epoch
position(s::State) = s.pos
velocity(s::State) = s.vel
state(s::State) = (s.pos, s.vel)
keplerian(s::State) = keplerian(position(s), velocity(s), grav_param(centralbody(s)))

struct KeplerianState{S,F,B,T,ET} <: AbstractState{S,F,B,T,ET}
    epoch::Epoch{S,ET}
    frame::F
    body::B
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    ν::T

    function KeplerianState(
        epoch::Epoch{S1,ET},
        a::A,
        e::E,
        i::AT,
        Ω::AT,
        ω::AT,
        ν::AT;
        scale::S2=timescale(epoch),
        frame::F=icrf,
        body::B=earth
    ) where {S1,S2,F,B,A,E,AT,ET}
        ep = Epoch{S2}(epoch)
        T = Union{A,E,AT}
        new{S2,F,B,T,ET}(ep, frame, body, a, e, i, Ω, ω, ν)
    end
end

function Base.getindex(s::KeplerianState, i::Int)
    i == 1 && return s.a
    i == 2 && return s.e
    i == 3 && return s.i
    i == 4 && return s.Ω
    i == 5 && return s.ω
    i == 6 && return s.ν
end

function (==)(s1::KeplerianState, s2::KeplerianState)
    return s1.epoch == s2.epoch &&
        s1.a == s2.a &&
        s1.e == s2.e &&
        s1.i == s2.i &&
        s1.Ω == s2.Ω &&
        s1.ω == s2.ω &&
        s1.ν == s2.ν
end

function Base.isapprox(
    s1::KeplerianState,
    s2::KeplerianState,
    atol::Real=0,
    rtol::Real=Base.rtoldefault(
        promote_leaf_eltypes(s1),
        promote_leaf_eltypes(s2),
        atol,
    ),
)
    return isapprox(s1.epoch, s2.epoch, atol=atol, rtol=rtol) &&
        isapprox(s1.a, s2.a, atol=atol, rtol=rtol) &&
        isapprox(s1.e, s2.e, atol=atol, rtol=rtol) &&
        isapprox(s1.i, s2.i, atol=atol, rtol=rtol) &&
        isapprox(s1.Ω, s2.Ω, atol=atol, rtol=rtol) &&
        isapprox(s1.ω, s2.ω, atol=atol, rtol=rtol) &&
        isapprox(s1.ν, s2.ν, atol=atol, rtol=rtol)
end

epoch(s::KeplerianState) = s.epoch
position(s::KeplerianState) = cartesian(s..., grav_param(body(s)))[1]
velocity(s::KeplerianState) = cartesian(s..., grav_param(body(s)))[2]
state(s::KeplerianState) = cartesian(s..., grav_param(body(s)))
keplerian(s::KeplerianState) = (s.a, s.e, s.i, s.Ω, s.ω, s.ν)

function KeplerianState(s::AbstractState)
    return KeplerianState(epoch(s), keplerian(s)...; frame=refframe(s), body=centralbody(s))
end
State(s::KeplerianState) = State(epoch(s), state(s)...; frame=refframe(s), body=body(s))

# struct ThreeBodyState{
#         F<:Frame,
#         T,
#         C1<:CelestialBody,
#         C2<:CelestialBody,
#     } <: AbstractState
#     ep::Epoch{T}
#     r::SVector{3,Float64}
#     v::SVector{3,Float64}
#
#     function ThreeBodyState(
#         ep::Epoch{T}, r, v, frame::Type{F}=GCRF,
#         primary::Type{C1}=Sun, secondary::Type{C2}=Earth
#     ) where {T,F<:Frame,C1<:CelestialBody,C2<:CelestialBody}
#         new{F, T, C1, C2}(ep, r, v)
#     end
# end
#
# frame(::ThreeBodyState{F}) where F<:Frame = F
# timescale(::ThreeBodyState{<:Frame, T}) where T = T
# primary(::ThreeBodyState{<:Frame,,C}) where C<:CelestialBody = C
# secondary(::ThreeBodyState{<:Frame,,<:CelestialBody,C}) where {
#     C<:CelestialBody} = C
