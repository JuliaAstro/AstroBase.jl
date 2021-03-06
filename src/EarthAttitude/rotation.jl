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

using ..Time: Epoch, SECONDS_PER_DAY, TDB, TT, UT1, julian_period
using ..Util: normalize2pi, sec_to_rad

export apparent_sidereal, earth_rotation_angle, equinoxes, mean_sidereal

function earth_rotation_angle(::IAU2000, ep::Epoch)
    t = julian_period(Float64, ep; scale=UT1)
    f = t % 1.0
    return mod2pi(2π * (f + 0.7790572732640 + 0.00273781191135448t))
end

include(joinpath("constants", "equinoxes.jl"))

function equinoxes(::IAU1994, ep::Epoch; scale=TDB)
    t = julian_period(Float64, ep; scale=scale, unit=centuries)

    ω₀ = sec_to_rad(@evalpoly(t, 450160.280, -482890.539, 7.455, 0.008))
    ω = mod2pi(ω₀ + (-5.0t % 1.0) * 2π)

    δψ, _ = nutation(iau1980, ep; scale=scale)
    ϵ₀ = obliquity(iau1980, ep; scale=scale)

    return δψ * cos(ϵ₀) + sec_to_rad(0.00264 * sin(ω) + 0.000063 * sin(ω + ω))
end

function equinoxes(::IAU2000, ep::Epoch; scale=TT)
    t = julian_period(Float64, ep; scale=scale, unit=centuries)

    l = fundamental(luna, t)
    lp = fundamental(sun, t)
    f = fundamental(luna, Longitude(), t)
    d = fundamental(luna, Elongation(), t)
    om = fundamental(luna, AscendingNode(), t)
    ve = fundamental(venus, t)
    le = fundamental(earth, t)
    pa = fundamental(t)

    s0 = 0.0
    s1 = 0.0

    for x in Iterators.reverse(E0)
        eq = (
            x.nl * l +
            x.nlp * lp +
            x.nf * f +
            x.nd * d +
            x.nom * om +
            x.nve * ve +
            x.nle * le +
            x.npa * pa
        )
        se, ce = sincos(eq)
        s0 += x.s * se + x.c * ce
    end

    for x in Iterators.reverse(E1)
        eq = (
            x.nl * l +
            x.nlp * lp +
            x.nf * f +
            x.nd * d +
            x.nom * om +
            x.nve * ve +
            x.nle * le +
            x.npa * pa
        )
        se, ce = sincos(eq)
        s1 += x.s * se + x.c * ce
    end

    return sec_to_rad(s0 + s1 * t)
end

function equinoxes(::IAU2000, ep::Epoch, ϵ, δψ; scale=TT)
    δψ * cos(ϵ) + equinoxes(iau2000, ep; scale=scale)
end

function equinoxes(iau::IAU2000Model, ep::Epoch; scale=TT)
    _, δϵ_pr = precession(iau2000, ep; scale=scale)
    ϵ = obliquity(iau1980, ep; scale=scale) + δϵ_pr
    δψ, _ = nutation(iau, ep; scale=scale)
    return equinoxes(iau2000, ep, ϵ, δψ; scale=scale)
end

function equinoxes(::IAU2006A, ep::Epoch)
    gst06a = apparent_sidereal(iau2006a, ep)
    gmst06 = mean_sidereal(iau2006, ep)
    return normalize2pi(gst06a - gmst06)
end

const SECS_TO_RAD = 7.272205216643039903848712e-5

function mean_sidereal(::IAU1982, ep::Epoch)
    t = julian_period(Float64, ep; scale=UT1, unit=centuries)
    f = julian_period(Float64, ep; scale=UT1) % 1.0 * SECONDS_PER_DAY
    gmst0 = @evalpoly(
        t,
        24110.54841 - SECONDS_PER_DAY / 2.0 + f,
        8640184.812866,
        0.093104,
        -6.2e-6,
    )
    return mod2pi(SECS_TO_RAD * gmst0)
end

function mean_sidereal(::IAU2000, ep::Epoch; scale=TT)
    t = julian_period(Float64, ep; scale=scale, unit=centuries)
    gmst0 = @evalpoly(
        t,
        0.014506,
        4612.15739966,
        1.39667721,
        -0.00009344,
        0.00001882,
    )
    return mod2pi(earth_rotation_angle(iau2000, ep) + sec_to_rad(gmst0))
end

function mean_sidereal(::IAU2006, ep::Epoch)
    t = julian_period(Float64, ep; scale=TT, unit=centuries)
    gmst0 = @evalpoly(
        t,
        0.014506,
        4612.156534,
        1.3915817,
        -0.00000044,
        -0.000029956,
        -0.0000000368,
    )
    return mod2pi(earth_rotation_angle(iau2000, ep) + sec_to_rad(gmst0))
end

function apparent_sidereal(::IAU1994, ep::Epoch)
    gmst82 = mean_sidereal(iau1982, ep)
    eqeq94 = equinoxes(iau1994, ep; scale=UT1)
    return mod2pi(gmst82 + eqeq94)
end

function apparent_sidereal(::IAU2000A, ep::Epoch)
    ee = equinoxes(iau2000a, ep)
    gmst = mean_sidereal(iau2000, ep)
    return mod2pi(gmst + ee)
end

function apparent_sidereal(::IAU2000B, ep::Epoch)
    ee = equinoxes(iau2000b, ep; scale=UT1)
    gmst = mean_sidereal(iau2000, ep; scale=UT1)
    return mod2pi(gmst + ee)
end

function apparent_sidereal(::IAU2006, ep::Epoch, rnpb)
    x, y = cip_coords(rnpb)
    s = cio_locator(iau2006, ep, x, y)
    era = earth_rotation_angle(iau2000, ep)
    eors = equation_of_origins(rnpb, s)
    return mod2pi(era - eors)
end

function apparent_sidereal(::IAU2006A, ep::Epoch)
    rnpb = precession_nutation_matrix(iau2006a, ep)
    return apparent_sidereal(iau2006, ep, rnpb)
end

