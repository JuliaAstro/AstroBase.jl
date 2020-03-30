#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
module TwoBody

using LinearAlgebra: norm, ×, ⋅
using StaticArrays: SVector
using ReferenceFrameRotations: angle_to_dcm

using Roots: find_zero, Order5, Order2

using ..Util: angle, azimuth, normalize_angle

export
    EccentricAnomaly,
    Elliptic,
    Hyperbolic,
    MeanAnomaly,
    Parabolic,
    TrueAnomaly,
    cartesian,
    eccentric_anomaly,
    elliptic,
    hyperbolic,
    isprograde,
    isretrograde,
    keplerian,
    mean_anomaly,
    parabolic,
    perifocal,
    semilatus,
    transform,
    true_anomaly

include("stumpff.jl")
include("kepler.jl")

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

iscircular(ecc, tol=1e-8) = isapprox(ecc, 0.0, atol=tol)
isequatorial(inc, tol=1e-8) = isapprox(abs(inc), 0.0, atol=tol)

function keplerian(pos, vel, µ, tol=1e-8)
    r = norm(pos)
    v = norm(vel)
    # Angular momentum
    h = pos × vel
    hm = norm(h)
    k = SVector(0.0, 0.0, 1.0)
    # Node vector
    node = k × h
    # Eccentricity vector
    e = ((v^2 - µ / r) * pos - (pos ⋅ vel) * vel) / µ
    ecc = norm(e)
    inc = angle(h, k)

    equatorial = isequatorial(inc, tol)
    circular = iscircular(ecc, tol)

    if circular
        # Semi-latus rectum
        a = hm^2 / µ
    else
        ξ = v^2 / 2 - µ / r
        a = -µ / (2ξ)
    end

    if equatorial && !circular
        Ω = 0.0
        # Longitude of pericenter
        ω = azimuth(e)
        ν = atan(h ⋅ (e × pos) / hm, pos ⋅ e)
    elseif !equatorial && circular
        Ω = azimuth(node)
        ω = 0.0
        # Argument of latitude
        ν = atan((pos ⋅ (h × node)) / hm, pos ⋅ node)
    elseif equatorial && circular
        Ω = 0.0
        ω = 0.0
        # True longitude
        ν = azimuth(pos)
    else
        if a > 0
            # Elliptic
            E_se = pos ⋅ vel / sqrt(μ * a)
            E_ce = r * v^2 / μ - 1
            ν = transform(elliptic, eccentric_anomaly, true_anomaly,
                          atan(E_se, E_ce), ecc)
        else
            # Hyperbolic
            E_sh = pos ⋅ vel / sqrt(-μ * a)
            E_ch = r * v^2 / μ - 1
            ν = transform(hyperbolic, eccentric_anomaly, true_anomaly,
                          log((E_ch + E_sh) / (E_ch - E_sh)) / 2, ecc)
        end
        Ω = azimuth(node)
        px = pos ⋅ node
        py = pos ⋅ (h × node) / hm
        ω = atan(py, px) - ν
    end
    Ω = mod2pi(Ω)
    ω = mod2pi(ω)
    ν = normalize_angle(ν, 0.0)

    a, ecc, inc, Ω, ω, ν
end

isprograde(inc) = inc < π/2 && !(inc ≈ π/2)
isretrograde(inc) = inc > π/2 && !(inc ≈ π/2)
ispolar(inc) = inc ≈ π/2

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
