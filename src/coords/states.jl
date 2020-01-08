using LinearAlgebra: promote_leaf_eltypes
import AstroTime: timescale

export State, KeplerianState, array, epoch, scale, frame, body, epoch_type, pos_type, vel_type

#################
# AbstractState #
#################

export AbstractState, epoch, frame

abstract type AbstractState{Scale, Frame, Body, T} <: AbstractArray{T, 1} end

function epoch end

timescale(s::AbstractState{Scale}) where {Scale} = Scale
frame(s::AbstractState{_S, Frame}) where {_S, Frame} = Frame
body(s::AbstractState{_S, _F, Body}) where {_S, _F, Body} = Body
array(s::AbstractState) = vcat(state(s)...)

Base.size(::AbstractState) = (6,)
Base.length(::AbstractState) = 6
Base.IndexStyle(::Type{<:AbstractState}) = IndexLinear()

function period(s::AbstractState)
    μ = grav_param(body(s))
    a, _ = keplerian(s)
    return period(a, μ)
end

struct State{Scale, Frame, Body, T, TP, TV} <: AbstractState{Scale, Frame, Body, T}
    epoch::Epoch{Scale, T}
    pos::SVector{3, TP}
    vel::SVector{3, TV}
    function State(epoch::Epoch{Scale, T}, pos, vel;
                   scale::TimeScale=timescale(epoch),
                   frame::Frame=icrf,
                   body::Body=earth) where {Scale, Frame, Body, T}
        TP = eltype(pos)
        TV = eltype(vel)
        ep = Epoch{scale}(epoch)
        new{scale::TimeScale, frame::AbstractFrame, body::CelestialBody, T, TP, TV}(ep, pos, vel)
    end
end

Base.getindex(s::State, i::Int) = getindex(s, Val(i))
Base.getindex(s::State, ::Val{1}) = s.pos[1]
Base.getindex(s::State, ::Val{2}) = s.pos[2]
Base.getindex(s::State, ::Val{3}) = s.pos[3]
Base.getindex(s::State, ::Val{4}) = s.vel[1]
Base.getindex(s::State, ::Val{5}) = s.vel[2]
Base.getindex(s::State, ::Val{6}) = s.vel[3]

function (==)(s1::State, s2::State)
    s1.epoch == s2.epoch && s1.pos == s2.pos && s1.vel == s2.vel
end

function Base.isapprox(s1::State, s2::State, atol::Real=0,
                       rtol::Real=Base.rtoldefault(promote_leaf_eltypes(s1), promote_leaf_eltypes(s2), atol))
    isapprox(s1.epoch, s2.epoch, atol=atol, rtol=rtol) &&
    isapprox(s1.pos, s2.pos, atol=atol, rtol=rtol) &&
    isapprox(s1.vel, s2.vel, atol=atol, rtol=rtol)
end

epoch_type(s::State{_S, _F, _B, T}) where {_S, _F, _B, T} = T
pos_type(s::State{_S, _F, _B, _T, TP}) where {_S, _F, _B, _T, TP} = TP
vel_type(s::State{_S, _F, _B, _T, _TP, TV}) where {_S, _F, _B, _T, _TP, TV} = TV

epoch(s::State) = s.epoch
position(s::State) = s.pos
velocity(s::State) = s.vel
state(s::State) = (s.pos, s.vel)
keplerian(s::State) = keplerian(position(s), velocity(s), grav_param(body(s)))

struct KeplerianState{Scale, Frame, Body, T} <: AbstractState{Scale, Frame, Body, T}
    epoch::Epoch{Scale, T}
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    ν::T
    function KeplerianState(epoch::Epoch{Scale, T}, a, e, i, Ω, ω, ν;
                            scale::TimeScale=timescale(epoch),
                            frame::Frame=icrf, body::Body=earth) where {Scale, Frame, Body, T}
        ep = Epoch{scale}(epoch)
        new{scale::TimeScale, frame::AbstractFrame, body::CelestialBody, T}(ep, a, e, i, Ω, ω, ν)
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

function Base.isapprox(s1::KeplerianState, s2::KeplerianState,
                       atol::Real=0,
                       rtol::Real=Base.rtoldefault(promote_leaf_eltypes(s1), promote_leaf_eltypes(s2), atol))
    isapprox(s1.epoch, s2.epoch, atol=atol, rtol=rtol) &&
    isapprox(s1.a, s2.a, atol=atol, rtol=rtol) &&
    isapprox(s1.e, s2.e, atol=atol, rtol=rtol) &&
    isapprox(s1.i, s2.i, atol=atol, rtol=rtol) &&
    isapprox(s1.Ω, s2.Ω, atol=atol, rtol=rtol) &&
    isapprox(s1.ω, s2.ω, atol=atol, rtol=rtol) &&
    isapprox(s1.ν, s2.ν, atol=atol, rtol=rtol)
end

Base.eltype(::KeplerianState{_S, _F, _B, T}) where {_S, _F, _B, T} = T

epoch(s::KeplerianState) = s.epoch
position(s::KeplerianState) = cartesian(s.a, s.e, s.i, s.Ω, s.ω, s.ν, grav_param(body(s)))[1]
velocity(s::KeplerianState) = cartesian(s.a, s.e, s.i, s.Ω, s.ω, s.ν, grav_param(body(s)))[2]
state(s::KeplerianState) = cartesian(s.a, s.e, s.i, s.Ω, s.ω, s.ν, grav_param(body(s)))
keplerian(s::KeplerianState) = (s.a, s.e, s.i, s.Ω, s.ω, s.ν)

KeplerianState(s::AbstractState) = KeplerianState(epoch(s), keplerian(s)...; frame=frame(s), body=body(s))
State(s::KeplerianState) = State(epoch(s), state(s)...; frame=frame(s), body=body(s))

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
