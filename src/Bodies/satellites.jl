#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

# Earth

@body luna 301 EarthSatellite parent=earth_barycenter _export=true

"""
    moon

`moon` is an alias for [`luna`](@ref).
"""
const moon = luna
export moon

# Mars

@body phobos 401 MarsSatellite parent=mars_barycenter _export=true
@body deimos 402 MarsSatellite parent=mars_barycenter _export=true

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
    name = Symbol(lowercase(body))
    id = 500 + i
    @eval begin
        @body $name $id JupiterSatellite parent=jupiter_barycenter _export=true
    end
end

@body dia 553 JupiterSatellite parent=jupiter_barycenter _export=true

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
    name = Symbol(lowercase(body))
    id = 600 + i
    @eval begin
        @body $name $id SaturnSatellite parent=saturn_barycenter _export=true
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
    name = Symbol(lowercase(body))
    id = 700 + i
    @eval begin
        @body $name $id UranusSatellite parent=uranus_barycenter _export=true
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
    name = Symbol(lowercase(body))
    id = 800 + i
    @eval begin
        @body $name $id NeptuneSatellite parent=neptune_barycenter _export=true
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
    name = Symbol(lowercase(body))
    id = 900 + i
    @eval begin
        @body $name $id PlutoSatellite parent=pluto_barycenter _export=true
    end
end

