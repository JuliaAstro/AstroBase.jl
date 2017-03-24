using AstroDynBase

import AstroDynBase: Rotation

export IAUFrame, body, Rotation, rotation

abstract type IAUFrame <: Frame end

for planet in PLANETS
    frame = Symbol("IAU", planet)
    @eval begin
        struct $frame <: IAUFrame end
        body(::Type{$frame}) = $planet
        export $frame
    end
end

function rotation(b::Type{C}, ep::TDBEpoch) where C<:CelestialBody
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
    m, δm = rotation(body(F), TDBEpoch(ep))
    M = zeros(6, 6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    M[4:6,1:3] = δm
    Rotation{GCRF,F}(M)
end

# IAU -> GCRF
function Rotation(::Type{F}, ::Type{GCRF}, ep::Epoch) where F<:IAUFrame
    m, δm = rotation(body(F), TDBEpoch(ep))
    M = zeros(6, 6)
    M[1:3,1:3] = m'
    M[4:6,4:6] = m'
    M[4:6,1:3] = δm'
    Rotation{F,GCRF}(M)
end
