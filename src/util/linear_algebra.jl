#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
export angle, azimuth, elevation, vector_azel

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

function azimuth(v)
    atan(v[2], v[1])
end

function elevation(v)
    asin(v[3] / norm(v))
end

function vector_azel(az, el)
    saz, caz = sincos(az)
    sel, cel = sincos(el)
    SVector(caz * cel, saz * cel, sel)
end
