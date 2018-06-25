module AstroBase

export sp00
const J000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const ARCSEC_TO_RAD = 4.848136811095359935899141e-6
# package code goes here
function sp00(jd1, jd2)
    t = (jd1 - J000 + jd2) / DAYS_PER_CENTURY
    sp = -47e-6 * t * ARCSEC_TO_RAD
    sp
end

end # module
