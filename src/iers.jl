struct CIRF <: Frame
struct TIRF <: Frame
struct ITRF <: Frame

function Rotation(::Type{CIRF}, ::Type{GCRF}, ep::Epoch)
    m = rotation_matrix(DATA.iau2000, TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m'
    M[4:6,4:6] = m'
    return M
end

function Rotation(::Type{GCRF}, ::Type{CIRF}, ep::Epoch)
    m = rotation_matrix(DATA.iau2000, TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    return M
end

function Rotation(::Type{TIRF}, ::Type{CIRF}, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = eraEra00(ut1.jd, ut1.jd1)
    rate = rotation_rate(EARTH, TDBEpoch(ep))
    m = rotation_matrix(3, era)
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    M[4:6,1:3] = rate_matrix(3, era, rate)
    return M
end

function Rotation(::Type{CIRF}, ::Type{TIRF}, ep::Epoch)
    ut1 = UT1Epoch(ep)
    era = eraEra00(ut1.jd, ut1.jd1)
    rate = rotation_rate(EARTH, TDBEpoch(ep))
    m = rotation_matrix(3, -era)
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    M[4:6,1:3] = rate_matrix(3, -era, -rate)
    return M
end

function polarmotion(data::IERSData, ep::TTEpoch)
    xp, yp = interpolate(data, ep)
    reshape(eraPom00(xp, yp, eraSp00(ep.jd, ep.jd1)), (3,3))
end

function precession_nutation(data::IERSData, ep::TTEpoch)
    dx, dy = interpolate(data, ep)
    x, y = eraXy06(ep.jd, ep.jd1)
    s = eraS06(ep.jd, ep.jd1, x, y)
    x += dx
    y += dy
    reshape(eraC2ixys(x, y, s), (3,3))
end

function Rotation(::Type{ITRF}, ::Type{TIRF}, ep::Epoch)
    m = rotation_matrix(DATA.polarmotion, TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m'
    M[4:6,4:6] = m'
    return M
end

function Rotation(::Type{TIRF}, ::Type{ITRF}, ep::Epoch)
    m = rotation_matrix(DATA.polarmotion, TTEpoch(ep))
    M = zeros(6,6)
    M[1:3,1:3] = m
    M[4:6,4:6] = m
    return M
end
