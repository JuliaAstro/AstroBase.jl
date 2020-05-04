#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module Bodies

import Base: parent

using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, add_vertex!, items

using ..Time: SECONDS_PER_DAY, SECONDS_PER_CENTURY, seconds, julian_period

export CelestialBody,
    NAIFId,
    SolarSystemBarycenter,
    Sun,
    along_orbit_radius,
    declination,
    declination_rate,
    ellipsoid,
    equatorial_radius,
    euler_angles,
    euler_rates,
    from_naifid,
    grav_param,
    mean_radius,
    naifid,
    parent,
    polar_radius,
    right_ascension,
    right_ascension_rate,
    rotation_angle,
    rotation_rate,
    ssb,
    subplanetary_radius,
    sun

include("abstract.jl")

const BODIES = ItemGraph{NAIFId, NAIFId}(SimpleGraph())


Base.show(io::IO, body::CelestialBody) = print(io, string(nameof(typeof(body))))

from_naifid(id::NAIFId) = from_naifid(Val(id))

path_ids(from::CelestialBody, to::CelestialBody) = items(BODIES, naifid(from), naifid(to))


struct SolarSystemBarycenter <: Barycenter end
const ssb = SolarSystemBarycenter()
Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
naifid(::SolarSystemBarycenter) = 0
from_naifid(::Val{0}) = ssb
add_vertex!(BODIES, 0)

struct Sun <: CelestialBody end
const sun = Sun()
parent(::Sun) = ssb
naifid(::Sun) = 10
from_naifid(::Val{10}) = sun
add_edge!(BODIES, 0, 10)


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
