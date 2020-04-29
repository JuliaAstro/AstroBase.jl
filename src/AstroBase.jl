#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module AstroBase

using Reexport

struct AstroException <: Exception
    msg::String
end

Base.showerror(io::IO, err::AstroException) = print(io, err.msg)

include(joinpath("Time", "Time.jl"))
include(joinpath("constants", "Constants.jl"))
include(joinpath("interfaces", "Interfaces.jl"))
include(joinpath("util", "Util.jl"))
include(joinpath("bodies", "Bodies.jl"))
include(joinpath("earth_attitude", "EarthAttitude.jl"))
include(joinpath("two_body", "TwoBody.jl"))
include(joinpath("ephemerides", "Ephemerides.jl"))
include(joinpath("frames", "Frames.jl"))
include(joinpath("coords", "Coords.jl"))

# include(joinpath("n_body", "NBody.jl"))
# include(joinpath("astrometry", "Astrometry.jl"))

@reexport using .Time
@reexport using .Constants
@reexport using .Interfaces
@reexport using .Util
@reexport using .Bodies
@reexport using .EarthAttitude
@reexport using .TwoBody
@reexport using .Ephemerides
@reexport using .Frames
@reexport using .Coords

# @reexport using .Astrometry

end # module
