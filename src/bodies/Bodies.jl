module Bodies

import Base: parent

using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, add_vertex!, items

export CelestialBody,
    NAIFId,
    SolarSystemBarycenter,
    Sun,
    from_naifid,
    grav_param,
    naifid,
    parent,
    ssb,
    sun

const NAIFId = Int
const bodies = ItemGraph{NAIFId, NAIFId}(SimpleGraph())

abstract type CelestialBody end

function grav_param end
function mean_radius end
function polar_radius end
function equatorial_radius end

Base.show(io::IO, body::CelestialBody) = print(io, string(nameof(typeof(body))))

from_naifid(id::NAIFId) = from_naifid(Val(id))

path_ids(from::CelestialBody, to::CelestialBody) = items(bodies, naifid(from), naifid(to))

abstract type Barycenter <: CelestialBody end

struct SolarSystemBarycenter <: Barycenter end
const ssb = SolarSystemBarycenter()
Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
naifid(::SolarSystemBarycenter) = 0
from_naifid(::Val{0}) = ssb

struct Sun <: CelestialBody end
const sun = Sun()
parent(::Sun) = ssb
naifid(::Sun) = 10
from_naifid(::Val{10}) = sun

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

include(joinpath(@__DIR__, "..", "..", "gen", "gm.jl"))

function __init__()
    # Sun and SSB
    add_vertex!(bodies, 0)
    add_edge!(bodies, 0, 10)

    # Planets and barycenters
    for i in 1:length(PLANET_NAMES)
        id = 100i + 99
        add_edge!(bodies, 0, i)
        add_edge!(bodies, i, id)
    end

    # Earth satellite
    add_edge!(bodies, 3, 301)

    # Mars satellites
    add_edge!(bodies, 4, 401)
    add_edge!(bodies, 4, 402)

    # Jupiter satellites
    for i in 1:length(JUPITER_SATELLITE_NAMES[1:end-1])
        id = 500 + i
        add_edge!(bodies, 5, id)
    end
    add_edge!(bodies, 5, 553)

    # Saturn satellites
    for i in 1:length(SATURN_SATELLITE_NAMES)
        id = 600 + i
        add_edge!(bodies, 6, id)
    end

    # Uranus satellites
    for i in 1:length(URANUS_SATELLITE_NAMES)
        id = 700 + i
        add_edge!(bodies, 7, id)
    end

    # Neptune satellites
    for i in 1:length(NEPTUNE_SATELLITE_NAMES)
        id = 800 + i
        add_edge!(bodies, 8, id)
    end

    # Pluto satellites
    for i in 1:length(PLUTO_SATELLITE_NAMES)
        id = 900 + i
        add_edge!(bodies, 9, id)
    end

    # Minor bodies
    add_edge!(bodies, 0, 9)
    add_edge!(bodies, 9, 999)
    add_edge!(bodies, 0, 2000001)
    add_edge!(bodies, 0, 2000002)
    add_edge!(bodies, 0, 2000004)
    add_edge!(bodies, 0, 2000016)
    add_edge!(bodies, 0, 2000021)
    add_edge!(bodies, 0, 2431010)
    add_edge!(bodies, 0, 2000433)
    add_edge!(bodies, 0, 2000511)
    add_edge!(bodies, 0, 9511010)
    add_edge!(bodies, 0, 2002867)
    add_edge!(bodies, 0, 2025143)
    add_edge!(bodies, 0, 1000093)
    add_edge!(bodies, 0, 1000005)
end

end
