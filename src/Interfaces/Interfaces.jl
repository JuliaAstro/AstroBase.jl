#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module Interfaces

using ..Time: SECONDS_PER_CENTURY, SECONDS_PER_DAY, julian_period, seconds

include("bodies.jl")
include("frames.jl")
include("ephemerides.jl")
include("coords.jl")

end

