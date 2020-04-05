#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

abstract type NaturalSatellite <: CelestialBody end

export deimos, Deimos, phobos, Phobos, luna, Luna, moon, charon, Charon

# Earth

struct Luna <: NaturalSatellite end
const luna = Luna()
const moon = luna
parent(::Luna) = earth_barycenter
naifid(::Luna) = 301
from_naifid(::Val{301}) = luna
add_edge!(BODIES, 3, 301)

# Mars

struct Phobos <: NaturalSatellite end
const phobos = Phobos()
parent(::Phobos) = mars_barycenter
naifid(::Phobos) = 401
from_naifid(::Val{401}) = phobos
add_edge!(BODIES, 4, 401)

struct Deimos <: NaturalSatellite end
const deimos = Deimos()
parent(::Deimos) = mars_barycenter
naifid(::Deimos) = 402
from_naifid(::Val{402}) = deimos
add_edge!(BODIES, 4, 402)

# Jupiter

const JUPITER_SATELLITE_NAMES = (
    "Io",
    "Europa",
    "Ganymede",
    "Callisto",
    "Amalthea",
    "Himalia",
    "Elara",
    "Pasiphae",
    "Sinope",
    "Lysithea",
    "Carme",
    "Ananke",
    "Leda",
    "Thebe",
    "Adrastea",
    "Metis",
    "Callirrhoe",
    "Themisto",
    "Magaclite",
    "Taygete",
    "Chaldene",
    "Harpalyke",
    "Kalyke",
    "Iocaste",
    "Erinome",
    "Isonoe",
    "Praxidike",
    "Autonoe",
    "Thyone",
    "Hermippe",
    "Aitne",
    "Eurydome",
    "Euanthe",
    "Euporie",
    "Orthosie",
    "Sponde",
    "Kale",
    "Pasithee",
    "Hegemone",
    "Mneme",
    "Aoede",
    "Thelxinoe",
    "Arche",
    "Kallichore",
    "Helike",
    "Carpo",
    "Eukelade",
    "Cyllene",
    "Kore",
    "Herse",
    "Dia",
)

for (i, body) in enumerate(JUPITER_SATELLITE_NAMES[1:end-1])
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 500 + i
    add_edge!(BODIES, 5, id)
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = jupiter_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end

struct Dia <: NaturalSatellite end
const dia = Dia()
parent(::Dia) = jupiter_barycenter
naifid(::Dia) = 553
from_naifid(::Val{553}) = dia
add_edge!(BODIES, 5, 553)
export dia, Dia

# Saturn

const SATURN_SATELLITE_NAMES = (
    "Mimas",
    "Enceladus",
    "Tethys",
    "Dione",
    "Rhea",
    "Titan",
    "Hyperion",
    "Iapetus",
    "Phoebe",
    "Janus",
    "Epimetheus",
    "Helene",
    "Telesto",
    "Calypso",
    "Atlas",
    "Prometheus",
    "Pandora",
    "Pan",
    "Ymir",
    "Paaliaq",
    "Tarvos",
    "Ijiraq",
    "Suttungr",
    "Kiviuq",
    "Mundilfari",
    "Albiorix",
    "Skathi",
    "Erriapus",
    "Siarnaq",
    "Thrymr",
    "Narvi",
    "Methone",
    "Pallene",
    "Polydeuces",
    "Daphnis",
    "Aegir",
    "Bebhionn",
    "Bergelmir",
    "Bestla",
    "Farbauti",
    "Fenrir",
    "Fornjot",
    "Hati",
    "Hyrrokkin",
    "Kari",
    "Loge",
    "Skoll",
    "Surtur",
    "Anthe",
    "Jarnsaxa",
    "Greip",
    "Tarqeq",
    "Aegaeon",
)

for (i, body) in enumerate(SATURN_SATELLITE_NAMES)
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 600 + i
    add_edge!(BODIES, 6, id)
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = saturn_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end

# Uranus

const URANUS_SATELLITE_NAMES = (
    "Ariel",
    "Umbriel",
    "Titania",
    "Oberon",
    "Miranda",
    "Cordelia",
    "Ophelia",
    "Bianca",
    "Cressida",
    "Desdemona",
    "Juliet",
    "Portia",
    "Rosalind",
    "Belinda",
    "Puck",
    "Caliban",
    "Sycorax",
    "Prospero",
    "Setebos",
    "Stephano",
    "Trinculo",
    "Francisco",
    "Margaret",
    "Ferdinand",
    "Perdita",
    "Mab",
    "Cupid",
)

for (i, body) in enumerate(URANUS_SATELLITE_NAMES)
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 700 + i
    add_edge!(BODIES, 7, id)
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = uranus_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end

# Neptune

const NEPTUNE_SATELLITE_NAMES = (
    "Triton",
    "Nereid",
    "Naiad",
    "Thalassa",
    "Despina",
    "Galatea",
    "Larissa",
    "Proteus",
    "Halimede",
    "Psamathe",
    "Sao",
    "Laomedeia",
    "Neso",
)

for (i, body) in enumerate(NEPTUNE_SATELLITE_NAMES)
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 800 + i
    add_edge!(BODIES, 8, id)
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = neptune_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end

# Pluto

const PLUTO_SATELLITE_NAMES = (
    "Charon",
    "Nix",
    "Hydra",
    "Kerberos",
    "Styx",
)

for (i, body) in enumerate(PLUTO_SATELLITE_NAMES)
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 900 + i
    add_edge!(BODIES, 9, id)
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = pluto_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end
