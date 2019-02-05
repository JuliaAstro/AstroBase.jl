module TwoBody

using LinearAlgebra: norm, ×, ⋅
using StaticArrays: SVector
using ReferenceFrameRotations: angle_to_dcm

using Roots: find_zero, Order5, Order2

import ..angle,
    ..azimuth,
    ..normalize_angle,
    ..vector_azel

export
    cartesian,
    keplerian,
    perifocal,
    semilatus,
    Elliptic,
    elliptic,
    Hyperbolic,
    hyperbolic,
    Parabolic,
    parabolic,
    MeanAnomaly,
    mean_anomaly,
    TrueAnomaly,
    true_anomaly,
    EccentricAnomaly,
    eccentric_anomaly,
    transform

semilatus(a, ecc) = ecc ≈ 0.0 ? a : a * (1.0 - ecc^2)

function cartesian(a, ecc, i, Ω, ω, ν, μ)
    p = semilatus(a, ecc)
    r_pqw, v_pqw = perifocal(p, ecc, ν, μ)
    M = angle_to_dcm(-ω, -i, -Ω,:ZXZ)

    M * r_pqw, M * v_pqw
end

function perifocal(p, ecc, ν, μ)
    sinν, cosν = sincos(ν)
    sqrt_μop = sqrt(μ/p)

    r_pqw = [cosν, sinν, zero(p)] .* p ./ (1 + ecc * cosν)
    v_pqw = [-sinν, ecc + cosν, zero(sqrt_μop)] .* sqrt_μop

    r_pqw, v_pqw
end

#= function cartesian(a, e, i, Ω, ω, v, μ) =#
#=     sinΩ, cosΩ = sincos(Ω) =#
#=     sinω, cosω = sincos(ω) =#
#=     sini, cosi = sincos(i) =#
#=  =#
#=     crcp = cosΩ * cosω =#
#=     crsp = cosΩ * sinω =#
#=     srcp = sinΩ * cosω =#
#=     srsp = sinΩ * sinω =#
#=  =#
#=     p = SVector(crcp - cosi * srsp, srcp + cosi * crsp, sini * sinω) =#
#=     q = SVector(-crsp - cosi * srcp, -srsp + cosi * crcp, sini * cosω) =#
#=  =#
#=     if a > 0 =#
#=         uME2 = (1 - e) * (1 + e) =#
#=         s1Me2 = sqrt(uME2) =#
#=         E = transform(elliptic, true_anomaly, eccentric_anomaly, v, e) =#
#=         cose = cos(E) =#
#=         sine = sin(E) =#
#=  =#
#=         x = a * (cose - e) =#
#=         y = a * sine * s1Me2 =#
#=         factor = sqrt(μ / a) / (1 - e * cose) =#
#=         ẋ = -sine * factor =#
#=         ẏ = cose * s1Me2 * factor =#
#=     else =#
#=         sinv = sin(v) =#
#=         cosv = cos(v) =#
#=         f = a * (1 - e^2) =#
#=         pos_factor = f / (1 + e * cosv) =#
#=         vel_factor = sqrt(μ / f) =#
#=  =#
#=         x = pos_factor * cosv =#
#=         y = pos_factor * sinv =#
#=         ẋ = -vel_factor * sinv =#
#=         ẏ = vel_factor * (e + cosv) =#
#=     end =#
#=  =#
#=     pos = x .* p .+ y .* q =#
#=     vel = ẋ .* p .+ ẏ .* q =#
#=  =#
#=     pos, vel =#
#= end =#

function keplerian(pos, vel, μ)
    h = pos × vel
    h2 = norm(h)^2
    k = SVector(0, 0, 1)
    i = angle(h, k)
    Ω = azimuth(k × h)

    r = norm(pos)
    r2 = r^2
    v2 = norm(vel)^2
    rv2_over_μ = r * v2 / μ

    a = r / (2 - rv2_over_μ)
    μa = μ * a

    if a > 0
        # Elliptic or circular
        e_se = pos ⋅ vel / sqrt(μa)
        e_ce = rv2_over_μ - 1
        ecc = sqrt(e_se^2 + e_ce^2)
        ν = transform(elliptic, eccentric_anomaly, true_anomaly,
                      atan(e_se, e_ce), ecc)
    else
        # Hyperbolic
        e_sh = pos ⋅ vel / sqrt(-μa)
        e_ch = rv2_over_μ - 1
        ecc = sqrt(1 - h2 / μa)
        ν = transform(hyperbolic, eccentric, true_anomaly,
                      log((e_ch + e_sh) / (e_ch - e_sh)) / 2, ecc)
    end

    node = vector_azel(Ω, 0.0)
    px = pos ⋅ node
    py = pos ⋅ (h × node) / sqrt(h2)
    ω = atan(py, px) - ν

    Ω = normalize_angle(Ω, π)
    ω = normalize_angle(ω, 0.0)
    ν = normalize_angle(ν, 0.0)

    a, ecc, i, Ω, ω, ν
end

abstract type OrbitType end

struct Elliptic <: OrbitType end
struct Hyperbolic <: OrbitType end
struct Parabolic <: OrbitType end

const elliptic = Elliptic()
const hyperbolic = Hyperbolic()
const parabolic = Parabolic()

abstract type Anomaly end

struct MeanAnomaly <: Anomaly end
struct TrueAnomaly <: Anomaly end
struct EccentricAnomaly <: Anomaly end

const mean_anomaly = MeanAnomaly()
const true_anomaly = TrueAnomaly()
const eccentric_anomaly = EccentricAnomaly()

kepler(::Elliptic, E, M, ecc) = E - ecc * sin(E) - M
kepler(::Hyperbolic, E, M, ecc) = -E + ecc * sinh(E) - M
kepler(::Parabolic, E, M, ecc) = mean_parabolic(E, ecc) - M

function mean_parabolic(E, ecc, tolerance=1e-16)
    x = (ecc - 1.0) / (ecc + 1.0) * E^2
    small = false
    S = 0.0
    k = 0
    while !small
        term = (ecc - 1.0 / (2.0k + 2.0)) * x^k
        small = abs(term) < tolerance
        S += term
        k += 1
    end
    sqrt(2.0 / (1.0 + ecc)) * E + sqrt(2.0 / (1.0 + ecc)^3) * E^3 * S
end

"""
    transform(from, to, a, ecc)

Transform anomaly `a` `from` one anomaly type `to` another for an orbit with eccentricity `ecc`.

### Arguments ###

- `from`, `to`: Anomaly types
    - `true_anomaly`
    - `eccentric_anomaly`
    - `mean_anomaly`
- `a`: Current value of the anomaly
- `ecc`: Eccentricity of the orbit

### Output ###

Returns the transformed anomaly.

### References ###

- Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
"""
function transform(from, to, a, ecc)
    ecc ≈ 1.0 && return transform(parabolic, from, to, a, 1.0)
    ecc > 1.0 && return transform(hyperbolic, from, to, a, ecc)

    transform(elliptic, from, to, a, ecc)
end

transform(::Parabolic, ::EccentricAnomaly, ::TrueAnomaly, a, _) = 2.0 * atan(a)
transform(::Parabolic, ::TrueAnomaly, ::EccentricAnomaly, a, _) = tan(a / 2.0)

function transform(::Elliptic, ::TrueAnomaly, ::EccentricAnomaly, ν, ecc)
    β = ecc / (1 + sqrt(1 - ecc^2))
    ν - 2 * atan(β * sin(ν) / (1 + β * cos(ν)))
end
function transform(::Elliptic, ::EccentricAnomaly, ::TrueAnomaly, E, ecc)
    β = ecc / (1 + sqrt((1 - ecc) * (1 + ecc)))
    E + 2 * atan(β * sin(E) / (1 - β * cos(E)))
end

function transform(::Hyperbolic, ::TrueAnomaly, ::EccentricAnomaly, a, ecc)
    log((sqrt(ecc + 1) + sqrt(ecc - 1) * tan(a / 2)) /
        (sqrt(ecc + 1) - sqrt(ecc - 1) * tan(a / 2)))
end
function transform(::Hyperbolic, ::EccentricAnomaly, ::TrueAnomaly, a, ecc)
    2 * atan((exp(a) * sqrt(ecc + 1) - sqrt(ecc + 1)) /
             (exp(a) * sqrt(ecc - 1) + sqrt(ecc - 1)))
end

function transform(::Elliptic, ::MeanAnomaly, ::EccentricAnomaly, M, ecc)
    find_zero(E -> kepler(elliptic, E, M, ecc), M, Order5())
end
transform(::Elliptic, ::EccentricAnomaly, ::MeanAnomaly, E, ecc) = kepler(elliptic, E, 0.0, ecc)

function transform(::Hyperbolic, ::MeanAnomaly, ::EccentricAnomaly, M, ecc)
    find_zero(E -> kepler(hyperbolic, E, M, ecc), asinh(M / ecc), Order2())
end
transform(::Hyperbolic, ::EccentricAnomaly, ::MeanAnomaly, E, ecc) = kepler(hyperbolic, E, 0.0, ecc)

function transform(::Parabolic, ::MeanAnomaly, ::EccentricAnomaly, M, ecc)
    B = 3.0 * M / 2.0
    A = (B + sqrt(1.0 + B^2))^(2.0 / 3.0)
    guess = 2 * A * B / (1 + A + A^2)
    find_zero(E -> kepler(parabolic, E, M, ecc), guess, Order5())
end
transform(::Parabolic, ::EccentricAnomaly, ::MeanAnomaly, E, ecc) = kepler(parabolic, E, 0.0, ecc)

function transform(::OrbitType, ::MeanAnomaly, ::TrueAnomaly, M, ecc, delta=1e-2)
    if ecc > 1 + delta
        typ = hyperbolic
    elseif ecc < 1 - delta
        typ = elliptic
    else
        typ = parabolic
    end

    E = transform(typ, mean_anomaly, eccentric_anomaly, M, ecc)
    transform(typ, eccentric_anomaly, true_anomaly, E, ecc)
end
function transform(::OrbitType, ::TrueAnomaly, ::MeanAnomaly, ν, ecc, delta=1e-2)
    if ecc > 1 + delta
        typ = hyperbolic
    elseif ecc < 1 - delta
        typ = elliptic
    else
        typ = parabolic
    end

    E = transform(typ, true_anomaly, eccentric_anomaly, ν, ecc)
    transform(typ, eccentric_anomaly, mean_anomaly, E, ecc)
end

end
