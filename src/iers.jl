using ERFA
import EarthOrientation: polarmotion, precession_nutation00
import AstronomicalTime: fjd1, fjd2
import AstroDynBase: Rotation

export CIRF, TIRF, ITRF, Rotation

struct CIRF <: Frame end
struct TIRF <: Frame end
struct ITRF <: Frame end

function Rotation(::Type{GCRF}, ::Type{CIRF}, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m'
    M[4:6,4:6] = m'
    Rotation{GCRF,CIRF}(M)
end

function Rotation(::Type{CIRF}, ::Type{GCRF}, ep::Epoch)
    m = precession_nutation(TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    Rotation{CIRF,GCRF}(M)
end

function Rotation(::Type{CIRF}, ::Type{TIRF}, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = eraEra00(fjd1(ut1), fjd2(ut1))
    rate = rotation_rate(Earth, TDBEpoch(ep))
    m = rotation_matrix(3, era)
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    M[4:6,1:3] = rate_matrix(3, era, rate)
    Rotation{CIRF,TIRF}(M)
end

function Rotation(::Type{TIRF}, ::Type{CIRF}, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = eraEra00(fjd1(ut1), fjd2(ut1))
    rate = rotation_rate(Earth, TDBEpoch(ep))
    m = rotation_matrix(3, -era)
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    M[4:6,1:3] = rate_matrix(3, -era, -rate)
    Rotation{TIRF,CIRF}(M)
end

function polarmotion(ep::TTEpoch)
    xp, yp = polarmotion(julian(ep))
    xp = dms2rad(0, 0, xp)
    yp = dms2rad(0, 0, yp)
    reshape(eraPom00(xp, yp, eraSp00(fjd1(ep), fjd2(ep))), (3,3))
end

function precession_nutation(ep::TTEpoch)
    dx, dy = precession_nutation00(julian(ep))
    x, y = eraXy06(fjd1(ep), fjd2(ep))
    s = eraS06(fjd1(ep), fjd2(ep), x, y)
    x += dms2rad(0, 0, dx/1000.0)
    y += dms2rad(0, 0, dy/1000.0)
    reshape(eraC2ixys(x, y, s), (3,3))
end

function Rotation(::Type{TIRF}, ::Type{ITRF}, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m'
    M[4:6,4:6] = m'
    Rotation{TIRF,ITRF}(M)
end

function Rotation(::Type{ITRF}, ::Type{TIRF}, ep::Epoch)
    m = polarmotion(TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    Rotation{ITRF,TIRF}(M)
end
