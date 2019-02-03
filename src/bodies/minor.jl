export MinorBody, pluto, Pluto, pluto_barycenter, PlutoBarycenter

abstract type MinorBody <: CelestialBody end

const MINOR_BODY_NAMES = (
    "Ceres",
    "Pallas",
    "Vesta",
    "Lutetia",
    "Ida",
    "Eros",
    "Davida",
    "Gaspra",
    "Steins",
    "Itokawa",
    "Tempel1",
    "Borrelly",
)

for body in MINOR_BODY_NAMES
    str_sym = Symbol(body)
    con_sym = Symbol(lowercase(body))
    @eval begin
        struct $str_sym <: MinorBody end
        const $con_sym = $str_sym()
        parent(::$str_sym) = ssb
        export $str_sym, $con_sym
    end
end

struct PlutoBarycenter <: Barycenter end
const pluto_barycenter = PlutoBarycenter()
Base.show(io::IO, ::PlutoBarycenter) = print(io, "Pluto Barycenter")
parent(::PlutoBarycenter) = ssb
naifid(::PlutoBarycenter) = 9
from_naifid(::Val{9}) = pluto_barycenter
add_edge!(bodies, 0, 9)

struct Pluto <: MinorBody end
const pluto = Pluto()
parent(::Pluto) = pluto_barycenter
naifid(::Pluto) = 999
from_naifid(::Val{999}) = pluto
add_edge!(bodies, 9, 999)

naifid(::Ceres) = 2000001
from_naifid(::Val{2000001}) = ceres
add_edge!(bodies, 0, 2000001)

naifid(::Pallas) = 2000002
from_naifid(::Val{2000002}) = pallas
add_edge!(bodies, 0, 2000002)

naifid(::Vesta) = 2000003
from_naifid(::Val{2000003}) = vesta
add_edge!(bodies, 0, 2000003)

naifid(::Lutetia) = 2000021
from_naifid(::Val{2000021}) = lutetia
add_edge!(bodies, 0, 2000021)

naifid(::Ida) = 2431010
from_naifid(::Val{2431010}) = ida
add_edge!(bodies, 0, 2431010)

naifid(::Eros) = 2000433
from_naifid(::Val{2000433}) = eros
add_edge!(bodies, 0, 2000433)

naifid(::Davida) = 2000511
from_naifid(::Val{2000511}) = davida
add_edge!(bodies, 0, 2000511)

naifid(::Gaspra) = 9511010
from_naifid(::Val{9511010}) = gaspra
add_edge!(bodies, 0, 9511010)

naifid(::Steins) = 2002867
from_naifid(::Val{2002867}) = steins
add_edge!(bodies, 0, 2002867)

naifid(::Itokawa) = 2025143
from_naifid(::Val{2025143}) = itokawa
add_edge!(bodies, 0, 2025143)

naifid(::Tempel1) = 1000093
from_naifid(::Val{1000093}) = tempel1
add_edge!(bodies, 0, 1000093)

naifid(::Borrelly) = 1000005
from_naifid(::Val{1000005}) = borrelly
add_edge!(bodies, 0, 1000005)

