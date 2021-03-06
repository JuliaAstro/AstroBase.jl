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

"""
    angle(v1, v2)

Returns the angle between vectors `v1` and `v2` in rad.

# Example

```jldoctest
julia> rad_to_deg(angle([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]))
90.0
```
"""
function Base.angle(v1, v2)
    normprod = norm(v1) * norm(v2)
    if normprod == 0.0
        throw(DomainError())
    else
        v1v2 = v1 ⋅ v2
        threshold = normprod * 0.9999
        if v1v2 >= -threshold && v1v2 <= threshold
            return acos(v1v2 / normprod)
        else
            v3n = norm(v1 × v2)
            return v1v2 >= 0.0 ? (v3n / normprod) : π - asin(v3n / normprod)
        end
    end
end

azimuth(v) = atan(v[2], v[1])

elevation(v) = asin(v[3] / norm(v))

function vector_azel(az, el)
    saz, caz = sincos(az)
    sel, cel = sincos(el)
    SVector(caz * cel, saz * cel, sel)
end

function angular_velocity(ψ, δψ, θ, δθ, ϕ, δϕ)
    Ω₁ = δψ * sin(θ) * sin(ϕ) - δθ * cos(ϕ)
    Ω₂ = δψ * sin(θ) * cos(ϕ) + δθ * sin(ϕ)
    Ω₃ = δψ * cos(θ) + δϕ
    [Ω₁, Ω₂, Ω₃]
end

@inbounds function frame(x)
    length(x) < 3 && throw(ArgumentError("Input vector must have three elements."))
    iszero(x) && return [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]

    x̂ = normalize(x[1:3])
    ŷ = Array{Float64}(undef, 3)
    ẑ = Array{Float64}(undef, 3)

    x2 = x̂ .^ 2
    _, idx = findmin(x2)
    s1, s2, s3 = circshift(1:3, 1-idx)
    f = sqrt(x2[s2] + x2[s3])

    ŷ[s1] = 0.0
    ŷ[s2] = -x̂[s3] / f
    ŷ[s3] = x̂[s2] / f
    ẑ[s1] = f
    ẑ[s2] = -x̂[s1] * ŷ[s3]
    ẑ[s3] = x̂[s1] * ŷ[s2]

    return x̂, ŷ, ẑ
end

function orthogonal(v)
    threshold = 0.6 * norm(v)
    iszero(threshold) && throw(ArgumentError("Norm of `v` is zero."))

    x, y, z = v
    if abs(x) <= threshold
        inverse = 1 / sqrt(y^2 + z^2)
        return [0.0, inverse * z, -inverse * y]
    elseif abs(y) <= threshold
        inverse = 1 / sqrt(x^2 + z^2)
        return [-inverse * z, 0.0, inverse * x]
    end
    inverse = 1 / sqrt(x^2 + y^2)
    return [inverse * y, -inverse * x, 0.0]
end

function plane_section(ellipsoid, plane_point, plane_normal)
    u = orthogonal(plane_normal)
    v = normalize(plane_normal × u)
    ue = u ./ ellipsoid
    ve = v ./ellipsoid
    pe = plane_point ./ ellipsoid
    alpha = ue ⋅ ue
    beta = ve ⋅ ve
    gamma = dot(ue, ve)
    delta = dot(pe, ue)
    epsilon = dot(pe, ve)
    zeta = dot([pe; 1.0], [pe; -1.0])
    # TODO: Check whether actually we need a more accurate dot product.
    # The deps of AccurateArithmetics.jl are not stable enough.
    # gamma = dot_oro(ue, ve)
    # delta = dot_oro(pe, ue)
    # epsilon = dot_oro(pe, ve)
    # zeta = dot_oro([pe; 1.0], [pe; -1.0])

    if abs(gamma) < floatmin(Float64)
        tan_theta = 0.0
    else
        bma = beta - alpha
        tan_theta = bma >= 0.0 ? -2gamma / (bma + sqrt(bma^2 + 4gamma^2)) : -2gamma / (bma - sqrt(bma^2 + 4gamma^2))
    end
    tan2 = tan_theta^2
    cos2 = 1 / (1 + tan2)
    sin2 = tan2 * cos2
    cos_sin = tan_theta * cos2
    co = sqrt(cos2)
    si = tan_theta * co

    # TODO: See above
    # denom = dot_oro([gamma, -alpha], [gamma, beta])
    # tau_c = dot_oro([beta, -gamma], [delta, epsilon]) / denom
    # nu_c = dot_oro([alpha, -gamma], [epsilon, delta]) / denom
    denom = dot([gamma, -alpha], [gamma, beta])
    tau_c = dot([beta, -gamma], [delta, epsilon]) / denom
    nu_c = dot([alpha, -gamma], [epsilon, delta]) / denom

    twogcs = 2 * gamma * cos_sin
    A = alpha * cos2 + beta * sin2 + twogcs
    B = alpha * sin2 + beta * cos2 - twogcs
    F = (alpha * tau_c + 2 * (gamma * nu_c + delta)) * tau_c + (beta * nu_c + 2 * epsilon) * nu_c + zeta
    l = sqrt(-F / A)
    m = sqrt(-F / B)

    isnan(l + m) && return nothing

    if l > m
        return plane_point .+ tau_c .* u .+ nu_c .* v, co .* u .+ si .* v, -si .* u .+ co .* v, l, m
    else
        return plane_point .+ tau_c .* u .+ nu_c .* v, si .* u .- co .* v, co .* u .+ si .* v, m, l
    end
end

"""
    distance(v1, v2)

Return the Euclidean (L2) distance between vectors `v1` and `v2`.

# Example

```jldoctest
julia> distance([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
1.4142135623730951
```
"""
distance(v1, v2) = norm(v1 .- v2)

function isinside(ellipsoid, point)
    scaled = ellipsoid .- point
    return scaled ⋅ scaled <= 1.0
end

function point_on_limb(ellipsoid, observer, outside)
    isinside(ellipsoid, observer) && throw(ArgumentError("`observer` is inside `ellipsoid`."))

    normal = observer × outside
    center, u, v, a, b = plane_section(ellipsoid, zeros(3), normal)
    a2 = a^2
    b2 = b^2
    delta = observer - center
    xo = delta ⋅ u
    yo = delta ⋅ v
    xo2 = xo^2
    yo2 = yo^2
    alpha = a2 * yo2 + b2 * xo2
    beta = a2 * b2 * xo
    gamma = a2 * a2 * (b2 - yo2)
    sq = sqrt(beta^2 - alpha * gamma)
    if beta > 0
        s = beta + sq
        xt1 = s / alpha
        xt2 = gamma / s
    else
        s = beta - sq
        xt1 = gamma / s
        xt2 = s / alpha
    end
    t1 = center .+ xt1 .* u .+ (b2 * (a2 - xt1 * xo) / (a2 * yo)) .* v
    t2 = center .+ xt2 .* u .+ (b2 * (a2 - xt2 * xo) / (a2 * yo)) .* v
    return distance(t1, outside) <= distance(t2, outside) ? t1 : t2
end


"""
    spherical_to_cartesian(theta, phi)

Convert spherical coordinates to Cartesian.

### Arguments ###

- `theta`: longitude angle in radians
- `phi`: latitude angle in radians

### Returns ###
- `x`: magnitude of projection on x axis
- `y`: magnitude of projection on y axis
- `z`: magnitude of projection on z axis

### References ###

- [ERFA](https://github.com/liberfa/erfa/blob/master/src/s2c.c)
"""
function spherical_to_cartesian(theta, phi)
    cp = cos(phi);
    x = cos(theta) * cp;
    y = sin(theta) * cp;
    z = sin(phi);

    return x, y, z
end
