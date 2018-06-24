module AstroBase

export a00, anp
const J2000 = 2451545.0
anp = x -> mod(x, 2pi) < 0? mod(x, 2pi)+2pi:mod(x, 2pi)

function a00(jd1, jd2)
    if jd1 < jd2
        d1 = jd1
        d2 = jd2
    else
        d1 = jd2
        d2 = jd1
    end
     t = d1 + d2 - J2000

     f = mod(d1, 1.0) + mod(d2, 1.0)

     anp(2pi * (f + 0.7790572732640 + 0.00273781191135448 * t))
end

end
