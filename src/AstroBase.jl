module AstroBase

export earth_rotation_angle
const J2000 = 2451545.0

"""
    earth_rotation_angle(jd1, jd2)

Return Earth rotation angle (radians) for a given UT1 2-part Julian Date (jd1, jd2).

# Example
```jldoctest
julia> earth_rotation_angle(2.4578265e6, 0.30434616919175345)
4.912208135094597
```
"""
function earth_rotation_angle(jd1, jd2)
    if jd1 < jd2
        d1 = jd1
        d2 = jd2
    else
        d1 = jd2
        d2 = jd1
    end
    t = d1 + d2 - J2000
    f = mod(d1, 1.0) + mod(d2, 1.0)
    mod2pi(2pi * (f + 0.7790572732640 + 0.00273781191135448 * t))
end

end
