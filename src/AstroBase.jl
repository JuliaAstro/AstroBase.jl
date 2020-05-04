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
include(joinpath("Constants", "Constants.jl"))
include(joinpath("Interfaces", "Interfaces.jl"))
include(joinpath("Util", "Util.jl"))
include(joinpath("Bodies", "Bodies.jl"))
include(joinpath("EarthAttitude", "EarthAttitude.jl"))
include(joinpath("TwoBody", "TwoBody.jl"))
include(joinpath("Ephemerides", "Ephemerides.jl"))
include(joinpath("Frames", "Frames.jl"))
include(joinpath("Coords", "Coords.jl"))

# include(joinpath("NBody", "NBody.jl"))
# include(joinpath("Astrometry", "Astrometry.jl"))

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
