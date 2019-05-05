export MinorBody, pluto, Pluto, pluto_barycenter, PlutoBarycenter

abstract type MinorBody <: CelestialBody end

const MINOR_BODY_NAMES = (
    "Ceres",
    "Pallas",
    "Vesta",
    "Psyche",
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

struct Pluto <: MinorBody end
const pluto = Pluto()
parent(::Pluto) = pluto_barycenter
naifid(::Pluto) = 999
from_naifid(::Val{999}) = pluto

naifid(::Ceres) = 2000001
from_naifid(::Val{2000001}) = ceres

naifid(::Pallas) = 2000002
from_naifid(::Val{2000002}) = pallas

naifid(::Vesta) = 2000004
from_naifid(::Val{2000004}) = vesta

naifid(::Psyche) = 2000016
from_naifid(::Val{2000016}) = psyche

naifid(::Lutetia) = 2000021
from_naifid(::Val{2000021}) = lutetia

naifid(::Ida) = 2431010
from_naifid(::Val{2431010}) = ida

naifid(::Eros) = 2000433
from_naifid(::Val{2000433}) = eros

naifid(::Davida) = 2000511
from_naifid(::Val{2000511}) = davida

naifid(::Gaspra) = 9511010
from_naifid(::Val{9511010}) = gaspra

naifid(::Steins) = 2002867
from_naifid(::Val{2002867}) = steins

naifid(::Itokawa) = 2025143
from_naifid(::Val{2025143}) = itokawa

naifid(::Tempel1) = 1000093
Base.show(io::IO, ::Tempel1) = print(io, "Tempel 1")
from_naifid(::Val{1000093}) = tempel1

naifid(::Borrelly) = 1000005
from_naifid(::Val{1000005}) = borrelly

