module Bodies

import Base: parent

using ItemGraphs: ItemGraph, add_edge!, add_vertex!

export CelestialBody, SolarSystemBarycenter, ssb, Sun, sun, parent, naifid, from_naifid

bodies = ItemGraph{Int}()

abstract type CelestialBody end
Base.show(io::IO, x::CelestialBody) = print(io, string(typeof(x)))

from_naifid(id::Int) = from_naifid(Val(id))

abstract type Barycenter <: CelestialBody end

struct ParentOfSSB <: Barycenter end

struct SolarSystemBarycenter <: Barycenter end
const ssb = SolarSystemBarycenter()
Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
naifid(::SolarSystemBarycenter) = 0
from_naifid(::Val{0}) = ssb
add_vertex!(bodies, 0)

struct Sun <: CelestialBody end
const sun = Sun()
parent(::Sun) = ssb
naifid(::Sun) = 10
from_naifid(::Val{10}) = sun
add_edge!(bodies, 0, 10)

include("planets.jl")
include("minor.jl")
include("satellites.jl")

end
