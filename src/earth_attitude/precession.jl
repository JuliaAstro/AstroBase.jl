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
using ..Util: sec_to_rad

export
    bias,
    bias_precession_matrix,
    fukushima_williams,
    fukushima_williams_matrix,
    precession

function bias(::IAU2000)
    δψ_b = sec_to_rad(-0.041775)
    δϵ_b = sec_to_rad(-0.0068192)
    δra₀ = sec_to_rad(-0.0146)
    return δψ_b, δϵ_b, δra₀
end

function precession(::IAU2000, ep::Epoch; scale=TT)
    t = julian_period(ep; scale=scale, unit=centuries, raw=true)
    precor = sec_to_rad(-0.29965)
    oblcor = sec_to_rad(-0.02524)
    return precor * t, oblcor * t
end

function fukushima_williams(::IAU2006, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)

    γ = @evalpoly(
        t,
        -0.052928,
        10.556378,
        0.4932044,
        -0.00031238,
        -0.000002788,
        0.0000000260,
    ) |> sec_to_rad
    ϕ = @evalpoly(
        t,
        84381.412819,
        -46.811016,
        0.0511268,
        0.00053289,
        -0.000000440,
        -0.0000000176,
    ) |> sec_to_rad
    ψ = @evalpoly(
        t,
        -0.041775,
        5038.481484,
        1.5584175,
        -0.00018522,
        -0.000026452,
        -0.0000000148,
    ) |> sec_to_rad
    ϵ = obliquity(iau2006, ep)
    return γ, ϕ, ψ, ϵ
end

function fukushima_williams_matrix(γ, ϕ, ψ, ϵ)
    return compose_rotation(
        angleaxis_to_dcm(γ, [0.0, 0.0, 1.0]),
        angle_to_dcm(ϕ, -ψ, -ϵ, :XZX),
    )
end

function bias_precession_matrix(::IAU2000, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)

    ϵ₀ = sec_to_rad(84381.448)

    δψ_b, δϵ_b, δra₀ = bias(iau2000)

    ψ_a77 = sec_to_rad(@evalpoly(t, 0.0, 5038.7784, -1.07259, -0.001147))
    ω_a77  = ϵ₀ + sec_to_rad(@evalpoly(t, 0.0, 0.0, 0.05127, -0.007726))
    χ_a   = sec_to_rad(@evalpoly(t, 0.0, 10.5526, -2.38064, -0.001125))

    δψ_pr, δϵ_pr = precession(iau2000, ep)
    ψ_a = ψ_a77 + δψ_pr
    ω_a  = ω_a77 + δϵ_pr

    rb = angle_to_dcm(δra₀, δψ_b * sin(ϵ₀), -δϵ_b, :ZYX)
    rp = compose_rotation(
        angleaxis_to_dcm(ϵ₀, [1.0, 0.0, 0.0]),
        angle_to_dcm(-ψ_a, -ω_a, χ_a, :ZXZ),
    )

    return rb, rp, compose_rotation(rb, rp)
end

