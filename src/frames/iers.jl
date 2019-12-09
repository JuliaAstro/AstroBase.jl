import EarthOrientation
using ReferenceFrameRotations: angleaxis_to_dcm, ddcm

import ..sec2rad
using AstroTime: value, julian, julian_twopart, TTEpoch, UT1Epoch, TDBEpoch
using ..EarthAttitude:
    celestial_to_intermediate,
    earth_rotation_angle,
    polar_motion,
    s06,
    tio_locator,
    xy06

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
    m = precession_nutation(TTEpoch(ep))
    Rotation{icrf, cirf}(m')
end

function Rotation(::CIRF, ::ICRF, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    Rotation{cirf, icrf}(m)
end

function Rotation(::CIRF, ::TIRF, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = earth_rotation_angle(value.(julian_twopart(ut1))...)
    rate = rotation_rate(earth, TDBEpoch(ep))
    m = angleaxis_to_dcm(era, [0, 0, 1])
    m′ = ddcm(m, [0, 0, rate])
    Rotation{cirf, tirf}(m, m′)
end

function Rotation(::TIRF, ::CIRF, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = earth_rotation_angle(value.(julian_twopart(ut1))...)
    rate = rotation_rate(earth, TDBEpoch(ep))
    m = angleaxis_to_dcm(-era, [0, 0, 1])
    m′ = ddcm(m, [0, 0, -rate])
    Rotation{tirf, cirf}(m, m′)
end

function polarmotion(ep::TTEpoch)
    xp, yp = EarthOrientation.polarmotion(value(julian(ep)))
    xp = sec2rad(xp)
    yp = sec2rad(yp)
    sp00 = tio_locator(value.(julian_twopart(ep))...)
    reshape(polar_motion(xp, yp, sp00), (3,3))
end

function precession_nutation(ep::TTEpoch)
    dx, dy = EarthOrientation.precession_nutation00(value(julian(ep)))
    jd1, jd2 = value.(julian_twopart(ep))
    x, y = xy06(jd1, jd2)
    s = s06(jd1, jd2, x, y)
    x += sec2rad(dx/1000.0)
    y += sec2rad(dy/1000.0)
    reshape(celestial_to_intermediate(x, y, s), (3,3))
end

function Rotation(::TIRF, ::ITRF, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    Rotation{tirf, itrf}(m')
end

function Rotation(::ITRF, ::TIRF, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    Rotation{itrf, tirf}(m)
end

