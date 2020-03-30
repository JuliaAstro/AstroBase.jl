#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# This file incorporates work covered by the following copyright and
# permission notice:
#
#   Copyright (C) 2013-2014, NumFOCUS Foundation.
#   All rights reserved.
#
#   This library is derived, with permission, from the International
#   Astronomical Union's "Standards of Fundamental Astronomy" library,
#   available from http://www.iausofa.org.
#
#   The ERFA version is intended to retain identical
#   functionality to the SOFA library, but made distinct through
#   different function and file names, as set out in the SOFA license
#   conditions. The SOFA original has a role as a reference standard
#   for the IAU and IERS, and consequently redistribution is permitted only
#   in its unaltered state. The ERFA version is not subject to this
#   restriction and therefore can be included in distributions which do not
#   support the concept of "read only" software.
#
#   Although the intent is to replicate the SOFA API (other than replacement of
#   prefix names) and results (with the exception of bugs; any that are
#   discovered will be fixed), SOFA is not responsible for any errors found
#   in this version of the library.
#
#   If you wish to acknowledge the SOFA heritage, please acknowledge that
#   you are using a library derived from SOFA, rather than SOFA itself.
#
#
#   TERMS AND CONDITIONS
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   1 Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2 Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#
#   3 Neither the name of the Standards Of Fundamental Astronomy Board, the
#      International Astronomical Union nor the names of its contributors
#      may be used to endorse or promote products derived from this software
#      without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

