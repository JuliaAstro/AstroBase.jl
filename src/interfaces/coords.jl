import AstroTime: timescale

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

abstract type AbstractCoord{S, F, B, T, N} <: AbstractArray{T, N} end
abstract type AbstractState{S, F, B, T} <: AbstractCoord{S, F, B, T, 1} end
abstract type AbstractTrajectory{S, F, B, T} <: AbstractCoord{S, F, B, T, 2} end

function epoch end

centralbody(::AbstractCoord{S, F, B}) where {S, F, B} = B
refframe(::AbstractCoord{S, F}) where {S, F} = F
timescale(::AbstractCoord{S}) where {S} = S
Base.eltype(::AbstractCoord{S, F, B, T}) where {S, F, B, T} = T

Base.size(::AbstractState) = (6,)
Base.length(::AbstractState) = 6
Base.IndexStyle(::Type{<:AbstractState}) = IndexLinear()

