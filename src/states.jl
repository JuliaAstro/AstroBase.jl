using AstronomicalTime

export AbstractState, State, ThreeBodyState

abstract type AbstractState end

struct State{T<:Timescale,F<:Frame,C<:CelestialBody} <: AbstractState
    ep::Epoch{T}
end

struct ThreeBodyState{T<:Timescale,F<:Frame,C1<:CelestialBody,C2<:CelestialBody} <: AbstractState
    ep::Epoch{T}
end
