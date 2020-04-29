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

using ReferenceFrameRotations: angle_to_dcm

using ..Time: Epoch, TT, centuries, julian_period
using ..Util: normalize2pi, sec_to_rad

export nutation, nutation_matrix

"""
    nutation(model, ep)

Return the nutation components for a given epoch and model.

# Arguments

- `model`: IAU model, one of: `iau1980`, `iau2000a`, `iau2000b`, `iau2006`
- `ep`: An epoch

# Output

- `δψ`: Nutation in longitude (radians)
- `δϵ`: Nutation in obliquity (radians)

# Example

```jldoctest
julia> ep = TTEpoch(2020, 1, 1)
2020-01-01T00:00:00.000 TT

julia> nutation(iau2006, ep)
(-7.996558232098883e-5, -8.25141288270117e-6)
```

# References

- [SOFA](http://www.iausofa.org/publications/sofa_pn.pdf)
"""
nutation

struct NutationLuniSolar
    # Coefficients of l,l',F,D,Om
    nl::Int
    nlp::Int
    nf::Int
    nd::Int
    nom::Int
    # longitude sin, t*sin, cos coefficients
    sp::Float64
    spt::Float64
    cp::Float64
    # obliquity cos, t*cos, sin coefficients
    ce::Float64
    cet::Float64
    se::Float64
end

struct NutationPlanetary
    # Coefficients of l, F, D and Omega
    nl::Int
    nf::Int
    nd::Int
    nom::Int
    # Coefficients of planetary longitudes
    nme::Int
    nve::Int
    nea::Int
    nma::Int
    nju::Int
    nsa::Int
    nur::Int
    nne::Int
    # Coefficient of general precession
    npa::Int
    # Longitude sin, cos coefficients
    sp::Int
    cp::Int
    # Obliquity sin, cos coefficients
    se::Int
    ce::Int
end

include(joinpath("constants", "nutation_1980.jl"))
include(joinpath("constants", "nutation_2000a.jl"))
include(joinpath("constants", "nutation_2000b.jl"))

function nutation(::IAU1980, ep; scale=TT)
    t = julian_period(Float64, ep; scale=scale, unit=centuries)

    el_poly = @evalpoly t 485866.733 715922.633 31.310 0.064
    el = normalize2pi(sec_to_rad(el_poly) + (1325.0 * t % 1.0) * 2π)
    elp_poly = @evalpoly t 1287099.804 1292581.224 -0.577 -0.012
    elp = normalize2pi(sec_to_rad(elp_poly) + (99.0 * t % 1.0) * 2π)
    f_poly = @evalpoly t 335778.877 295263.137 -13.257 0.011
    f = normalize2pi(sec_to_rad(f_poly) + (1342.0 * t % 1.0) * 2π)
    d_poly = @evalpoly t 1072261.307 1105601.328 -6.891 0.019
    d = normalize2pi(sec_to_rad(d_poly) + (1236.0 * t % 1.0) * 2π)
    om_poly = @evalpoly t 450160.280 -482890.539 7.455 0.008
    om = normalize2pi(sec_to_rad(om_poly) + (-5.0 * t % 1.0) * 2π)

    dp = 0.0
    de = 0.0

    for x in Iterators.reverse(NUTATION_1980)
        arg = x.nl * el + x.nlp * elp + x.nf  * f + x.nd  * d + x.nom * om
        sarg, carg = sincos(arg)

        s = x.sp + x.spt * t
        c = x.ce + x.cet * t
        dp += s * sarg
        de += c * carg
    end

    δψ = sec_to_rad(dp * 1e-4)
    δϵ = sec_to_rad(de * 1e-4)

    return δψ, δϵ
end

function nutation(::IAU2000A, ep::Epoch; scale=TT)
    t = julian_period(Float64, ep; scale=scale, unit=centuries)

    el = fundamental(luna, t)
    f  = fundamental(luna, Longitude(), t)
    om = fundamental(luna, AscendingNode(), t)

    elp_as = @evalpoly t 1287104.79305 129596581.0481 -0.5532 0.000136 -0.00001149
    elp = sec_to_rad(elp_as % ARCSECONDS_IN_CIRCLE)
    d_as = @evalpoly t 1072260.70369 1602961601.2090 -6.3706 0.006593 -0.00003169
    d   = sec_to_rad(d_as % ARCSECONDS_IN_CIRCLE)

    dpls = 0.0
    dels = 0.0

    for x in Iterators.reverse(NUTATION_2000A_LS)
        arg = mod2pi(x.nl * el + x.nlp * elp + x.nf * f + x.nd * d + x.nom * om)
        sarg, carg = sincos(arg)

        dpls += (x.sp + x.spt * t) * sarg + x.cp * carg
        dels += (x.ce + x.cet * t) * carg + x.se * sarg
    end

    al   = mod2pi(2.35555598 + 8328.6914269554 * t)
    af   = mod2pi(1.627905234 + 8433.466158131 * t)
    ad   = mod2pi(5.198466741 + 7771.3771468121 * t)
    aom  = mod2pi(2.18243920 - 33.757045 * t)
    alne = mod2pi(5.321159000 + 3.8127774000 * t)

    apa  = fundamental(t)
    alme = fundamental(mercury, t)
    alve = fundamental(venus, t)
    alea = fundamental(earth, t)
    alma = fundamental(mars, t)
    alju = fundamental(jupiter, t)
    alsa = fundamental(saturn, t)
    alur = fundamental(uranus, t)

    dppl = 0.0
    depl = 0.0

    for x in Iterators.reverse(NUTATION_2000A_PL)
        arg = mod2pi(x.nl * al +
                     x.nf * af   +
                     x.nd * ad   +
                     x.nom * aom  +
                     x.nme * alme +
                     x.nve * alve +
                     x.nea * alea +
                     x.nma * alma +
                     x.nju * alju +
                     x.nsa * alsa +
                     x.nur * alur +
                     x.nne * alne +
                     x.npa * apa)
        sarg, carg = sincos(arg)

        dppl += x.sp * sarg + x.cp * carg
        depl += x.se * sarg + x.ce * carg
    end

    # Luni-solar terms
    δψ_ls = sec_to_rad(dpls * 1e-7)
    δϵ_ls = sec_to_rad(dels * 1e-7)

    # Planetary terms
    δψ_pl = sec_to_rad(dppl * 1e-7)
    δϵ_pl = sec_to_rad(depl * 1e-7)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end

function nutation(::IAU2000B, ep::Epoch; scale=TT)
    t = julian_period(Float64, ep; scale=scale, unit=centuries)

    el  = (485868.249036 + (1717915923.2178) * t) % ARCSECONDS_IN_CIRCLE |> sec_to_rad
    elp = (1287104.79305 + (129596581.0481) * t) % ARCSECONDS_IN_CIRCLE |> sec_to_rad
    f   = (335779.526232 + (1739527262.8478) * t) % ARCSECONDS_IN_CIRCLE |> sec_to_rad
    d   = (1072260.70369 + (1602961601.2090) * t) % ARCSECONDS_IN_CIRCLE |> sec_to_rad
    om  = (450160.398036 + (-6962890.5431) * t) % ARCSECONDS_IN_CIRCLE |> sec_to_rad

    dp = 0.0
    de = 0.0

    for x in Iterators.reverse(NUTATION_2000B)
        arg = mod2pi(x.nl * el + x.nlp * elp + x.nf * f + x.nd * d + x.nom * om)
        sarg, carg = sincos(arg)

        dp += (x.sp + x.spt * t) * sarg + x.cp * carg
        de += (x.ce + x.cet * t) * carg + x.se * sarg
    end

    # Luni-solar terms
    δψ_ls = sec_to_rad(dp * 1e-7)
    δϵ_ls = sec_to_rad(de * 1e-7)

    # Fixed offsets in lieu of planetary terms
    δψ_pl = sec_to_rad(-0.135 * 1e-3)
    δϵ_pl = sec_to_rad(0.388 * 1e-3)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end

function nutation(::IAU2006A, ep::Epoch)
    t = julian_period(Float64, ep; scale=TT, unit=centuries)

    fj2 = -2.7774e-6 * t
    δψ, δϵ = nutation(iau2000a, ep)

    δψ += δψ * (0.4697e-6 + fj2)
    δϵ += δϵ * fj2

    return δψ, δϵ
end

function nutation_matrix(ϵ, δψ, δϵ)
    angle_to_dcm(ϵ, -δψ, -(ϵ + δϵ), :XZX)
end

function nutation_matrix(iau, ep)
    δψ, δϵ = nutation(iau, ep)
    ϵ = obliquity(iau, ep)
    return nutation_matrix(ϵ, δψ, δϵ)
end

