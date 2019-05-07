import ERFA
import EarthOrientation: polarmotion, precession_nutation00
using ReferenceFrameRotations: angleaxis_to_dcm, ddcm

using AstroTime: value, julian, julian_twopart, TTEpoch, UT1Epoch, TDBEpoch
import ..sec2rad

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

struct TIRF <: RotatingFrame end
const tirf = TIRF()
from_sym(::Val{:TIRF}) = tirf

struct ITRF <: RotatingFrame end
const itrf = ITRF()
from_sym(::Val{:ITRF}) = itrf

function Rotation(::ICRF, ::CIRF, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    Rotation{ICRF, CIRF}(m')
end

function Rotation(::CIRF, ::ICRF, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    Rotation{CIRF, ICRF}(m)
end

function Rotation(::CIRF, ::TIRF, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = ERFA.era00(value.(julian_twopart(ut1))...)
    rate = rotation_rate(earth, TDBEpoch(ep))
    # m = rotation_matrix(3, era)
    m = angleaxis_to_dcm(era, [0, 0, 1])
    # δm = rate_matrix(3, era, rate)
    m′ = ddcm(m, [0, 0, rate])
    Rotation{CIRF, TIRF}(m, m′)
end

function Rotation(::TIRF, ::CIRF, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = ERFA.era00(value.(julian_twopart(ut1))...)
    rate = rotation_rate(earth, TDBEpoch(ep))
    # m = rotation_matrix(3, -era)
    m = angleaxis_to_dcm(-era, [0, 0, 1])
    # δm = rate_matrix(3, -era, -rate)
    m′ = ddcm(m, [0, 0, -rate])
    Rotation{TIRF, CIRF}(m, m′)
end

function polarmotion(ep::TTEpoch)
    xp, yp = polarmotion(value(julian(ep)))
    xp = sec2rad(xp)
    yp = sec2rad(yp)
    reshape(ERFA.pom00(xp, yp, ERFA.sp00(value.(julian_twopart(ep))...)), (3,3))
end

function precession_nutation(ep::TTEpoch)
    dx, dy = precession_nutation00(value(julian(ep)))
    jd1, jd2 = value.(julian_twopart(ep))
    x, y = ERFA.xy06(jd1, jd2)
    s = ERFA.s06(jd1, jd2, x, y)
    x += sec2rad(dx/1000.0)
    y += sec2rad(dy/1000.0)
    reshape(ERFA.c2ixys(x, y, s), (3,3))
end

function Rotation(::TIRF, ::ITRF, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    Rotation{TIRF,ITRF}(m')
end

function Rotation(::ITRF, ::TIRF, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    Rotation{ITRF,TIRF}(m)
end

