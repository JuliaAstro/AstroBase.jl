#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module Bodies

using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, add_vertex!, items
using ..Time: SECONDS_PER_DAY, SECONDS_PER_CENTURY, seconds, julian_period

import Base: parent

export CelestialBody, NAIFId, SolarSystemBarycenter, Sun,
    along_orbit_radius, declination, declination_rate, ellipsoid, equatorial_radius,
    euler_angles, euler_rates, from_naifid, grav_param, mean_radius, naifid,
    parent, polar_radius, right_ascension, right_ascension_rate,
    rotation_angle, rotation_rate, ssb, subplanetary_radius, sun,
    @body

include("abstract.jl")
include("graph.jl")
include("planets.jl")
include("minor.jl")
include("satellites.jl")

const ALL_NAMES = Tuple([["Sun", "SolarSystemBarycenter"]
                         collect(PLANET_NAMES .* "Barycenter");
                         collect(PLANET_NAMES);
                         ["Luna", "Phobos", "Deimos"];
                         collect(JUPITER_SATELLITE_NAMES);
                         collect(SATURN_SATELLITE_NAMES);
                         collect(URANUS_SATELLITE_NAMES);
                         collect(NEPTUNE_SATELLITE_NAMES);
                         collect(MINOR_BODY_NAMES);
                         ["Pluto", "PlutoBarycenter"];
                         collect(PLUTO_SATELLITE_NAMES)])

gm = joinpath(@__DIR__, "..", "..", "gen", "gm.jl")
isfile(gm) || error("`gm.jl` has not been generated, yet. Run `gen/gen_gm.jl`")
include(gm)

pck = joinpath(@__DIR__, "..", "..", "gen", "pck.jl")
isfile(pck) || error("`pck.jl` has not been generated, yet. Run `gen/gen_pck.jl`")
include(pck)

end
