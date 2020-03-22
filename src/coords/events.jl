using LinearAlgebra: norm

using ..Bodies: CelestialBody, sun, mean_radius, ellipsoid
using ..Ephemerides: AbstractEphemeris, position
using ..Frames: Rotation, icrf, iau_earth
using ..TwoBody: isretrograde
using ..Util: angle, point_on_limb

abstract type Detector end
abstract type DiscreteDetector <: Detector end
abstract type IntervalDetector <: Detector end

struct Event{S,TT,T<:Detector}
    type::Symbol
    epoch::Epoch{S,TT}
    detector::T
end

struct Apocenter <: DiscreteDetector end

function detect(::Apocenter, t, tra)
    y = tra(t)
    el = keplerian(y[1:3], y[4:6], grav_param(centralbody(tra)))
    ano = el[6]
    if ano > pi/2
        ano = abs(ano - pi)
    elseif ano < -pi/2
        ano = -abs(ano + pi)
    end
    return isretrograde(el[3]) ? ano : -ano
end

struct Pericenter <: DiscreteDetector end

function detect(::Pericenter, t, tra)
    y = tra(t)
    el = keplerian(y[1:3], y[4:6], grav_param(centralbody(tra)))
    return isretrograde(el[3]) ? -el[6] : el[6]
end

struct Eclipse{T<:AbstractEphemeris, C1<:CelestialBody, C2<:CelestialBody} <: IntervalDetector
    eph::T
    occulting::C1
    occulted::C2
    total::Bool
end

function detect(ecl::Eclipse, t, tra)
    find_eclipse(ecl, epoch(tra) + t, tra(t)[1:3], centralbody(tra))
end

function find_eclipse(ecl::Eclipse, ep, pos, center)
    # FIXME: Use proper rotation instead of hard-coded values
    rot = Rotation(icrf, iau_earth, ep)
    pocc = position(ecl.eph, ep, center, ecl.occulting)
    pted, _ = rot(position(ecl.eph, ep, ecl.occulting, ecl.occulted), zeros(3))
    psat, _ = rot(pos .- pocc, zeros(3))
    plimb = point_on_limb(ellipsoid(ecl.occulting), psat, pted)
    ps = psat .- pted
    pi = psat .- plimb
    ang = angle(ps, psat)
    rs = asin(mean_radius(ecl.occulted) / norm(ps))
    isnan(rs) && return π
    ro = angle(pi, psat)
    val = ecl.total ? (ang - ro + rs) : (ang - ro - rs)
    return -val
end

