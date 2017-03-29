using AstronomicalTime

import AstroDynBase: AbstractState, keplerian, velocity, Rotation, period
import Base: show, isapprox
import Base.Operators: ==

export State, ThreeBodyState, period
export timescale, frame, body, primary, secondary, keplerian, radius, velocity,
    epoch, isapprox, ==

struct State{
        F<:Frame,
        T<:Timescale,
        C<:CelestialBody,
    } <: AbstractState
    epoch::Epoch{T}
    r::VectorKM
    v::VectorKPS

    function State(ep::Epoch{T}, r, v,
        frame::Type{F}=GCRF, body::Type{C}=Earth) where {
        F<:Frame, T<:Timescale, C<:CelestialBody}
        new{F,T,C}(ep, r, v)
    end
end

function State(ep::Epoch{T},
    sma, ecc, inc, node, peri, ano,
    frame::Type{F}=GCRF, body::Type{C}=Earth) where {
    T<:Timescale,F<:Frame,C<:CelestialBody}
    r, v = cartesian(sma, ecc, inc, node, peri, ano, μ(body))
    State(ep, r, v, frame, body)
end

function State(ep::Epoch{T}, rv,
    frame::Type{F}=GCRF, body::Type{C}=Earth) where {
    T<:Timescale,F<:Frame,C<:CelestialBody}
    State(ep, rv[1], rv[2], frame, body)
end

radius(s::State) = s.r
velocity(s::State) = s.v
rv(s::State) = [s.r; s.v]
epoch(s::State) = s.epoch
keplerian(s::State) = keplerian(radius(s), velocity(s), μ(body(s)))
function period(s::State)
    ele = keplerian(s)
    period(ele[1], μ(body(s)))
end

semimajor(s::State) = keplerian(s)[1]
eccentricity(s::State) = keplerian(s)[2]
inclination(s::State) = keplerian(s)[3]
ascendingnode(s::State) = keplerian(s)[4]
argofpericenter(s::State) = keplerian(s)[5]
trueanomaly(s::State) = keplerian(s)[6]
frame(::State{F}) where F<:Frame = F
const _frame = frame
timescale(::State{<:Frame, T}) where T<:Timescale = T
const _timescale = timescale
body(::State{<:Frame, <:Timescale, C}) where C<:CelestialBody = C
const _body = body
(rot::Rotation)(s::State) = rot(radius(s), velocity(s))

function State(s::State{F1, T1, C1};
    frame::Type{F2}=_frame(s), timescale::Type{T2}=_timescale(s),
    body::Type{C2}=_body(s)) where {F1<:Frame, F2<:Frame,
    T1<:Timescale, T2<:Timescale, C1<:CelestialBody, C2<:CelestialBody}
    convert(State{F2, T2, C2}, s)
end

function (==)(s1::State{F, T, C}, s2::State{F, T, C}) where {
    F<:Frame, T<:Timescale, C<:CelestialBody}
    s1.epoch == s2.epoch && s1.r == s2.r && s1.v == s2.v
end
(==)(s1::State{<:Frame, <:Timescale, <:CelestialBody},
    s2::State{<:Frame, <:Timescale, <:CelestialBody}) = false

function isapprox(s1::State{F, T, C}, s2::State{F, T, C}) where {
    F<:Frame, T<:Timescale, C<:CelestialBody}
    s1.epoch ≈ s2.epoch && s1.r ≈ s2.r && s1.v ≈ s2.v
end
isapprox(s1::State{<:Frame, <:Timescale, <:CelestialBody},
    s2::State{<:Frame, <:Timescale, <:CelestialBody}) = false

function show(io::IO, s::State)
    sma, ecc, inc, node, peri, ano = keplerian(s)
    print(io, "State{",
    frame(s), ", ",
    timescale(s), ", ",
    body(s), "}\n",
    " Epoch: ", s.epoch, "\n",
    " Frame: ", frame(s), "\n",
    " Body:  ", body(s), "\n\n",
    " rx: ", s.r[1], "\n",
    " ry: ", s.r[2], "\n",
    " rz: ", s.r[3], "\n",
    " vx: ", s.v[1], "\n",
    " vy: ", s.v[2], "\n",
    " vz: ", s.v[3], "\n\n",
    " a: ", sma, "\n",
    " e: ", ecc, "\n",
    " i: ", °(inc), "\n",
    " ω: ", °(node), "\n",
    " Ω: ", °(peri), "\n",
    " ν: ", °(ano))
end

struct ThreeBodyState{
        F<:Frame,
        T<:Timescale,
        C1<:CelestialBody,
        C2<:CelestialBody,
    } <: AbstractState
    ep::Epoch{T}
    r::VectorKM
    v::VectorKPS

    function ThreeBodyState(
        ep::Epoch{T}, r, v, frame::Type{F}=GCRF,
        primary::Type{C1}=Sun, secondary::Type{C2}=Earth
    ) where {T<:Timescale,F<:Frame,C1<:CelestialBody,C2<:CelestialBody}
        new{F, T, C1, C2}(ep, r, v)
    end
end

frame(::ThreeBodyState{F}) where F<:Frame = F
timescale(::ThreeBodyState{<:Frame, T}) where T<:Timescale = T
primary(::ThreeBodyState{<:Frame,<:Timescale,C}) where C<:CelestialBody = C
secondary(::ThreeBodyState{<:Frame,<:Timescale,<:CelestialBody,C}) where {
    C<:CelestialBody} = C
