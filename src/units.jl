using Unitful

import Unitful: km, s, kg, °, rad,
    𝐋, 𝐓, Length, Time

export km, s, kg, °, rad, kps, dms2rad, rad2dms, VectorKM, VectorKPS

const kps = km/s

@derived_dimension Velocity 𝐋/𝐓

VectorKM = typeof([1.0]*km)
VectorKPS = typeof([1.0]*kps)

function dms2rad(deg, arcmin, arcsec)
    deg2rad(deg + arcmin/60 + arcsec/3600)
end

function rad2dms(rad)
    d = rad2deg(rad)
    deg = trunc(d)
    arcmin = trunc((d-deg)*60)
    arcsec = (d-deg-arcmin/60)*3600
    return deg, arcmin, arcsec
end
