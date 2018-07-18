module Bodies

import Base: parent

export CelestialBody, SolarSystemBarycenter, ssb, Sun, sun, parent, naifid

abstract type CelestialBody end

abstract type Barycenter <: CelestialBody end

struct ParentOfSSB <: Barycenter end

struct SolarSystemBarycenter <: Barycenter end
const ssb = SolarSystemBarycenter()
Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
parent(::SolarSystemBarycenter) = ParentOfSSB()
naifid(::SolarSystemBarycenter) = 0

struct Sun <: CelestialBody end
const sun = Sun()
Base.show(io::IO, ::Sun) = print(io, "Sun")
parent(::Sun) = ssb
naifid(::Sun) = 10

include("planets.jl")
include("minor.jl")
include("satellites.jl")

end
