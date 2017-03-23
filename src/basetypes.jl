export AbstractState
export Frame, GCRF

abstract type AbstractState end

abstract type Frame end

struct GCRF <: Frame end
