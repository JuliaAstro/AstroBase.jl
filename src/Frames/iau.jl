#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

using ..Bodies
using ..Time: TDBEpoch

using ReferenceFrameRotations: angle_to_dcm, ddcm

for name in Bodies.ALL_NAMES
    endswith(name, "Barycenter") && continue

    frame = Symbol("IAU", name)
    cname = Symbol("iau_", lowercase(name))

    @eval begin
        @frame $cname type=$frame parent=icrf rotating=true
        export $cname, $frame

        function Rotation(::ICRF, ::$frame, ep::Epoch)
            tdb = TDBEpoch(ep)
            body = $(Symbol(name))()
            angles = euler_angles(body, tdb)
            rates = euler_rates(body, tdb)
            m = angle_to_dcm(angles..., :ZXZ)
            return Rotation(icrf, $cname, m, rates)
        end

        Rotation(::$frame, ::ICRF, ep::Epoch) = inv(Rotation(icrf, $cname, ep))
    end
end

