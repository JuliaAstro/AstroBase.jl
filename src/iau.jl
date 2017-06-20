using AstroDynBase

import AstroDynBase: Rotation, isrotating

export IAUFrame, body, Rotation, rotation

abstract type IAUFrame <: Frame end

for body in [PLANETS; SATELLITES; MINOR_BODIES]
    frame = Symbol("IAU", body)
    @eval begin
        struct $frame <: IAUFrame end
        body(::Type{$frame}) = $body
        isrotating(::Type{$frame}) = true
        export $frame
    end
end

function rotation_matrices(b::Type{C}, ep::TDBEpoch) where C<:CelestialBody
    α = right_ascension(b, ep)
    δα = right_ascension_rate(b, ep)
    δ = declination(b, ep)
    δδ = declination_rate(b, ep)
    ω = rotation_angle(b, ep)
    δω = rotation_rate(b, ep)
    ϕ = α + π/2
    χ = π/2 - δ

    m = rotation_matrix(313, ϕ, χ, ω)
    δm = rate_matrix(313, ϕ, δα, χ, -δδ, ω, δω)
    return m, δm
end

# GCRF -> IAU
function Rotation(::Type{GCRF}, ::Type{F}, ep::Epoch) where F<:IAUFrame
    m, δm = rotation_matrices(body(F), TDBEpoch(ep))
    Rotation{GCRF,F}(m, δm)
end

# IAU -> GCRF
function Rotation(::Type{F}, ::Type{GCRF}, ep::Epoch) where F<:IAUFrame
    m, δm = rotation_matrices(body(F), TDBEpoch(ep))
    Rotation{F,GCRF}(m', δm')
end
