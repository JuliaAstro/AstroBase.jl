export AbstractState, AbstractTrajectory
export Frame, GCRF, epoch, isrotating

abstract type AbstractState end
abstract type AbstractTrajectory end

abstract type Frame end

isrotating(::Type{F}) where {F<:Frame} = false

struct GCRF <: Frame end

Base.show(io::IO, ::Type{F}) where F<:Frame = print(io, Base.datatype_name(F))

function epoch() end
