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
        export $con_sym
    end
end

struct PlutoBarycenter <: Barycenter end
const pluto_barycenter = PlutoBarycenter()
Base.show(io::IO, ::PlutoBarycenter) = print(io, "Pluto Barycenter")
parent(::PlutoBarycenter) = ssb
naifid(::PlutoBarycenter) = 901

struct Pluto <: MinorBody end
const pluto = Pluto()
parent(::Pluto) = pluto_barycenter
naifid(::Pluto) = 999

naifid(::Ceres) = 2000001
naifid(::Pallas) = 2000002
naifid(::Vesta) = 2000003
naifid(::Lutetia) = 2000021
naifid(::Ida) = 2431010
naifid(::Eros) = 2000433
naifid(::Davida) = 2000511
naifid(::Gaspra) = 9511010
naifid(::Steins) = 2002867
naifid(::Itokawa) = 2025143
naifid(::Tempel1) = 1000093
naifid(::Borrelly) = 1000005

