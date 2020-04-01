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
using ..Bodies
using ..Util: sec_to_rad

export fundamental

const ARCSECONDS_IN_CIRCLE = 1296000.0

function fundamental(::Sun, t)
    p = @evalpoly(t, 1287104.793048, 129596581.0481, -0.5532, 0.000136, -0.00001149)
    return sec_to_rad(p % ARCSECONDS_IN_CIRCLE)
end

function fundamental(::Luna, t)
    p = @evalpoly(t, 485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470)
    return sec_to_rad(p % ARCSECONDS_IN_CIRCLE)
end

struct Longitude end

function fundamental(::Luna, ::Longitude, t)
    p = @evalpoly(t, 335779.526232, 1739527262.8478, -12.7512, -0.001037, 0.00000417)
    return sec_to_rad(p % ARCSECONDS_IN_CIRCLE)
end

struct Elongation end

function fundamental(::Luna, ::Elongation, t)
    p = @evalpoly(t, 1072260.703692, 1602961601.2090, -6.3706, 0.006593, -0.00003169)
    return sec_to_rad(p % ARCSECONDS_IN_CIRCLE)
end

struct AscendingNode end

function fundamental(::Luna, ::AscendingNode, t)
    p = @evalpoly(t, 450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939)
    return sec_to_rad(p % ARCSECONDS_IN_CIRCLE)
end

fundamental(::Mercury, t) = mod2pi(4.402608842 + 2608.7903141574t)
fundamental(::Venus, t) = mod2pi(3.176146697 + 1021.3285546211t)
fundamental(::Earth, t) = mod2pi(1.753470314 + 628.3075849991t)
fundamental(::Mars, t) = mod2pi(6.203480913 + 334.0612426700t)
fundamental(::Jupiter, t) = mod2pi(0.599546497 + 52.9690962641t)
fundamental(::Saturn, t) = mod2pi(0.874016757 + 21.3299104960t)
fundamental(::Uranus, t) = mod2pi(5.481293872 + 7.4781598567t)
fundamental(::Neptune, t) = mod2pi(5.311886287 + 3.8133035638t)
fundamental(t) = @evalpoly(t, 0.0, 0.024381750, 0.00000538691)

