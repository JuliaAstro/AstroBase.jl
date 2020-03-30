#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
using ..Time: Epoch, TT, centuries, julian_period
using ..Util: sec2rad

export obliquity


"""
    obliquity(model, ep)

Return the mean obliquity of the ecliptic for a given epoch and model.

# Arguments

- `model`: IAU model, one of: `iau1980`, `iau2006`
- `ep`: An epoch

# Output

Returns the angle between the ecliptic and mean equator of date in radians.

# Example

```jldoctest
julia> ep = TTEpoch(2020, 1, 1)
2020-01-01T00:00:00.000 TT

julia> obliquity(iau2006, ep)
0.40904718953841473
```

# References

- [SOFA](http://www.iausofa.org/publications/sofa_pn.pdf)
"""
obliquity

function obliquity(::IAU1980Model, ep::Epoch; scale=TT)
    t = julian_period(ep; scale=scale, unit=centuries, raw=true)
    obl = @evalpoly(t, 84381.448, -46.8150, -0.00059, 0.001813)
    return sec2rad(obl)
end

function obliquity(::IAU2006Model, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)
    obl = @evalpoly(t, 84381.406, -46.836769, -0.0001831, 0.00200340,
                    -0.000000576, -0.0000000434)
    return sec2rad(obl)
end

