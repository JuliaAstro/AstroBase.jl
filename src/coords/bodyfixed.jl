#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
export RAC

using ..Frames: AbstractFrame, ICRF

struct RAC <: AbstractFrame end

function Rotation(::RAC, ::ICRF, r, v)
    normal = cross(r, v)
    normal = normal / norm(normal)
    tangential = v / norm(v)
    orthogonal = cross(v, normal)
    orthogonal = orthogonal / norm(orthogonal)
    m = hcat(orthogonal, tangential, normal)
    Rotation{RAC(), icrf}(m)
end
