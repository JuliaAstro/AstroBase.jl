using ERFA
import EarthOrientation: polarmotion, precession_nutation00
import AstroDynBase: Rotation

export CIRF, TIRF, ITRF, Rotation

struct CIRF <: Frame end
struct TIRF <: Frame end
struct ITRF <: Frame end

function Rotation(::Type{GCRF}, ::Type{CIRF}, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    Rotation{GCRF,CIRF}(m')
end

function Rotation(::Type{CIRF}, ::Type{GCRF}, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    Rotation{CIRF,GCRF}(m)
end

function Rotation(::Type{CIRF}, ::Type{TIRF}, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = eraEra00(julian1_strip(ut1), julian2_strip(ut1))
    rate = rotation_rate(Earth, TDBEpoch(ep))
    m = rotation_matrix(3, era)
    δm = rate_matrix(3, era, rate)
    Rotation{CIRF,TIRF}(m, δm)
end

function Rotation(::Type{TIRF}, ::Type{CIRF}, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = eraEra00(julian1_strip(ut1), julian2_strip(ut1))
    rate = rotation_rate(Earth, TDBEpoch(ep))
    m = rotation_matrix(3, -era)
    δm = rate_matrix(3, -era, -rate)
    Rotation{TIRF,CIRF}(m, δm)
end

function polarmotion(ep::TTEpoch)
    xp, yp = polarmotion(julian_strip(ep))
    xp = dms2rad(0, 0, xp)
    yp = dms2rad(0, 0, yp)
    reshape(eraPom00(xp, yp, eraSp00(julian1_strip(ep), julian2_strip(ep))), (3,3))
end

function precession_nutation(ep::TTEpoch)
    dx, dy = precession_nutation00(julian_strip(ep))
    x, y = eraXy06(julian1_strip(ep), julian2_strip(ep))
    s = eraS06(julian1_strip(ep), julian2_strip(ep), x, y)
    x += dms2rad(0, 0, dx/1000.0)
    y += dms2rad(0, 0, dy/1000.0)
    reshape(eraC2ixys(x, y, s), (3,3))
end

function Rotation(::Type{TIRF}, ::Type{ITRF}, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    Rotation{TIRF,ITRF}(m')
end

function Rotation(::Type{ITRF}, ::Type{TIRF}, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    Rotation{ITRF,TIRF}(m)
end
