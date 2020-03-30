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

