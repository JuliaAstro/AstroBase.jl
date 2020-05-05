const BODIES = ItemGraph{NAIFId, NAIFId}(SimpleGraph())

path_ids(from::CelestialBody, to::CelestialBody) = items(BODIES, naifid(from), naifid(to))

macro body()
end

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

