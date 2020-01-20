import Base: position

export AbstractEphemeris,
    position,
    position!,
    velocity,
    velocity!,
    state,
    state!

abstract type AbstractEphemeris end

function position! end
function velocity end
function velocity! end
function state end
function state! end

