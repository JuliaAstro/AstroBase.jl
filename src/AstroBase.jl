module AstroBase

export tio_locator
const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const ARCSEC_TO_RAD = 4.848136811095359935899141e-6
# package code goes here
"""
    tio_locator(jd1, jd2)

The TIO locator s', positioning the Terrestrial Intermediate Origin
on the equator of the Celestial Intermediate Pole.

'''
jldoctest
julia> AstroBase.tio_locator(2.4578265e6, 0.30434616919175345)
-3.9189245827947945e-11
'''
"""
function tio_locator(jd1, jd2)
    t = (jd1 - J2000 + jd2) / DAYS_PER_CENTURY
    -47e-6 * t * ARCSEC_TO_RAD
end

end # module
