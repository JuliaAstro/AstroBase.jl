module AstroBase

export tio_locator, sec2rad, rad2sec
const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
# package code goes here
"""
    sec2rad(sec::Real)

Returns radians for given seconds.

# Example

```jldoctest
julia> sec2rad(3600 * 30)
0.5235987755982988
```
"""
sec2rad(sec::Real) = deg2rad(sec / 3600)

"""
    rad2sec(rad::Real)
Returns angle (radians) for given seconds.

# Example

```jldoctest
julia> rad2sec(0.5235987755982988)
107999.99999999999
```
"""
rad2sec(rad::Real) = rad2deg(rad) * 3600

"""
    tio_locator(jd1, jd2)

Returns TIO locator s' position for a given 2-part Julian date (jd1, jd2).

# Example

'''jldoctest
julia> AstroBase.tio_locator(2.4578265e6, 0.30434616919175345)
-3.9189245827947945e-11
'''
"""
function tio_locator(jd1, jd2)
    t = (jd1 - J2000 + jd2) / DAYS_PER_CENTURY
    -47e-6 * t * sec2rad(1)
end

end # module
