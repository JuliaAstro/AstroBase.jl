using ..TwoBody: isretrograde

struct Event{Scale, TType}
    epoch::Epoch{Scale, TType}
    typ::Symbol
end

abstract type Detector end

struct Apocenter <: Detector end

function detect(::Apocenter, t, tra)
    y = tra(t)
    el = keplerian(y[1:3], y[4:6], grav_param(body(tra)))
    ano = el[6]
    if ano > pi/2
        ano = abs(ano - pi)
    elseif ano < -pi/2
        ano = -abs(ano + pi)
    end
    return isretrograde(el[3]) ? ano : -ano
end

struct Pericenter <: Detector end

function detect(::Pericenter, t, tra)
    y = tra(t)
    el = keplerian(y[1:3], y[4:6], grav_param(body(tra)))
    return isretrograde(el[3]) ? -el[6] : el[6]
end
