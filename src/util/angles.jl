#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
export
    angle,
    days_to_hms,
    deg_to_dms,
    deg_to_rad,
    deg_to_sec,
    dms_to_deg,
    dms_to_rad,
    hms_to_days,
    normalize2pi,
    rad_to_deg,
    rad_to_dms,
    rad_to_sec,
    sec_to_deg,
    sec_to_rad

const ARCSECONDS_PER_DEGREE = 3600

"""
    deg_to_rad(x)

Convert `x` from degrees to radians.

# Example

```jldoctest
julia> deg_to_rad(90)
1.5707963267948966
"""
deg_to_rad(x) = deg2rad(x)

"""
    rad_to_deg(x)

Convert `x` from radians to degrees.

# Example

```jldoctest
julia> rad_to_deg(π/2)
90
"""
rad_to_deg(x) = rad2deg(x)

"""
    sec_to_rad(x)

Convert `x` from arcseconds to radians.

# Example

```jldoctest
julia> sec_to_rad(324000)
1.5707963267948966
```
"""
sec_to_rad(x) = deg_to_rad(x / ARCSECONDS_PER_DEGREE)

"""
    rad_to_sec(x)

Convert `x` from radians to arcseconds.

# Example

```jldoctest
julia> rad_to_sec(π/2)
324000.0
```
"""
rad_to_sec(x) = rad_to_deg(x) * ARCSECONDS_PER_DEGREE

"""
    sec_to_deg(x)

Convert `x` from arcseconds to degrees.

# Example

```jldoctest
julia> sec_to_deg(324000)
90.0
```
"""
sec_to_deg(x) = x / ARCSECONDS_PER_DEGREE

"""
    deg_to_sec(x)

Convert `x` from degrees to arcseconds.

# Example

```jldoctest
julia> deg_to_sec(90)
324000.0
```
"""
deg_to_sec(x) = x * ARCSECONDS_PER_DEGREE

"""
    deg_to_dms(x)

Split `x` into degrees, arcminutes, and arcseconds.

# Example

```jldoctest
julia> deg_to_dms(90.314159)
(90, 18, 50.97240000001307)
```
"""
function deg_to_dms(x)
    d, res = divrem(x, 1.0)
    m, res = divrem(res * 60.0, 1.0)
    s = res * 60.0
    return Int(d), Int(m), s
end

"""
    dms_to_deg(d, m, s)

Convert degrees `d`, arcminutes `m`, and arcseconds `s` to degrees.

# Example

```jldoctest
julia> dms_to_deg(90, 18, 50.97240000001307)
90.314159
```
"""
dms_to_deg(d, m, s) = ((s / 60.0) + m) / 60.0 + d


"""
    rad_to_dms(x)

Convert `x` from radians into degrees, arcminutes, and arcseconds.

# Example

```jldoctest
julia> rad_to_dms(π/7)
(25, 42, 51.42857142857508)
```
"""
rad_to_dms(x) = rad_to_deg(x) |> deg_to_dms

"""
    dms_to_rad(d, m, s)

Convert degrees `d`, arcminutes `m`, and arcseconds `s` to radians.

# Example

```jldoctest
julia> dms_to_rad(25, 42, 51.42857142857508)
0.4487989505128276
```
"""
dms_to_rad(d, m, s) = dms_to_deg(d, m, s) |> deg_to_rad

"""
    days_to_hms(days)

Split `days` into hours, minutes, and seconds.

# Example

```jldoctest
julia> days_to_hms(0.314159)
(7, 32, 23.337600000000265)
```
"""
function days_to_hms(x)
    h, res = divrem(x * 24.0, 1.0)
    m, res = divrem(res * 60.0, 1.0)
    s = res * 60.0
    return Int(h), Int(m), s
end

"""
    hms_to_days(h, m, s)

Convert hours `h`, minutes `m`, and seconds `s` to days.

# Example

```jldoctest
julia> hms_to_days(7, 32, 50.2)
0.3144699074074074
```
"""
hms_to_days(h, m, s) = (((s / 60.0) + m) / 60.0 + h) / 24.0

"""
    normalize2pi(angle[, center=0.0])

Normalize `angle` to be in the interval `[center-π, center+π]`.
"""
function normalize2pi(angle, center=0.0)
    angle - 2π * floor((angle + π - center) / 2π)
end

