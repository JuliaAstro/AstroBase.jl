export AbstractState, AbstractTrajectory
export Frame, GCRF, epoch, isrotating

abstract type AbstractState end
abstract type AbstractTrajectory end

abstract type Frame end

isrotating(::Type{F}) where {F<:Frame} = false

struct GCRF <: Frame end

function epoch() end
