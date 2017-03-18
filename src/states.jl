using AstronomicalTime

export AbstractState, State, ThreeBodyState

abstract type AbstractState end

struct State{
        T<:Timescale,
        F<:Frame,
        C<:CelestialBody,
    } <: AbstractState
    ep::Epoch{T}
    r::VectorKM
    v::VectorKPS

    function State(ep::Epoch{T}, r, v, frame::Type{F}=GCRF, body::Type{C}=Earth) where {T<:Timescale,F<:Frame,C<:CelestialBody}
        new{T,F,C}(ep, r, v)
    end
end

struct ThreeBodyState{
        T<:Timescale,
        F<:Frame,
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
        new{T,F,C1,C2}(ep, r, v)
    end
end
