#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
import EarthOrientation
using ReferenceFrameRotations: angleaxis_to_dcm, ddcm

using ..EarthAttitude:
    celestial_to_intermediate,
    cio_locator,
    cip_coords,
    earth_rotation_angle,
    iau2000,
    iau2006,
    polar_motion,
    tio_locator
using ..Time: UTC, julian_period
using ..Util: sec2rad

export
    CIRF,
    ITRF,
    TIRF,
    cirf,
    itrf,
    tirf

struct CIRF <: InertialFrame end
const cirf = CIRF()
from_sym(::Val{:CIRF}) = cirf
add_edge!(FRAMES, :ICRF, :CIRF)

struct TIRF <: RotatingFrame end
const tirf = TIRF()
from_sym(::Val{:TIRF}) = tirf
add_edge!(FRAMES, :CIRF, :TIRF)

struct ITRF <: RotatingFrame end
const itrf = ITRF()
from_sym(::Val{:ITRF}) = itrf
add_edge!(FRAMES, :TIRF, :ITRF)

function Rotation(::ICRF, ::CIRF, ep::Epoch)
    m = precession_nutation(ep)
    return Rotation{icrf, cirf}(m)
end

function Rotation(::CIRF, ::ICRF, ep::Epoch)
    m = precession_nutation(ep)
    return Rotation{cirf, icrf}(m')
end

function Rotation(::CIRF, ::TIRF, ep::Epoch)
    era = earth_rotation_angle(iau2000, ep)
    rate = rotation_rate(earth, ep)
    m = angleaxis_to_dcm(era, [0, 0, 1])
    m′ = ddcm(m, [0, 0, rate])
    return Rotation{cirf, tirf}(m, m′)
end

function Rotation(::TIRF, ::CIRF, ep::Epoch)
    era = earth_rotation_angle(iau2000, ep)
    rate = rotation_rate(earth, ep)
    m = angleaxis_to_dcm(-era, [0, 0, 1])
    m′ = ddcm(m, [0, 0, -rate])
    return Rotation{tirf, cirf}(m, m′)
end

function polarmotion(ep::Epoch)
    jd = julian_period(ep; scale=UTC, origin=:julian, raw=true)
    xp, yp = EarthOrientation.polarmotion(jd)
    xp = sec2rad(xp)
    yp = sec2rad(yp)
    sp00 = tio_locator(iau2000, ep)
    return polar_motion(iau2000, xp, yp, sp00)
end

function precession_nutation(ep::Epoch)
    jd = julian_period(ep; scale=UTC, origin=:julian, raw=true)
    dx, dy = EarthOrientation.precession_nutation00(jd)
    x, y = cip_coords(iau2006, ep)
    s = cio_locator(iau2006, ep, x, y)
    x += sec2rad(dx/1000.0)
    y += sec2rad(dy/1000.0)
    return celestial_to_intermediate(x, y, s)
end

function Rotation(::TIRF, ::ITRF, ep::Epoch)
    m = polarmotion(ep)
    Rotation{tirf, itrf}(m)
end

function Rotation(::ITRF, ::TIRF, ep::Epoch)
    m = polarmotion(ep)
    Rotation{itrf, tirf}(m')
end

