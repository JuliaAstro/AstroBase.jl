abstract type NaturalSatellite <: CelestialBody end

export deimos, Deimos, phobos, Phobos, luna, Luna, moon, charon, Charon

# Mars

struct Phobos <: NaturalSatellite end
const phobos = Phobos()
Base.show(io::IO, ::Phobos) = print(io, "Phobos")
parent(::Phobos) = mars_barycenter
naifid(::Phobos) = 401

struct Deimos <: NaturalSatellite end
const deimos = Deimos()
Base.show(io::IO, ::Deimos) = print(io, "Deimos")
parent(::Deimos) = mars_barycenter
naifid(::Deimos) = 402

# Earth

struct Luna <: NaturalSatellite end
const luna = Luna()
const moon = luna
Base.show(io::IO, ::Luna) = print(io, "Luna")
parent(::Luna) = earth_barycenter
naifid(::Luna) = 301

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
    "Harpalke",
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

for (i, body) in enumerate(JUPITER_SATELLITE_NAMES)
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 500 + i
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        Base.show(io::IO, ::$typ) = print(io, $body)
        parent(::$typ) = jupiter
        export $sym, $typ
    end
end

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
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        Base.show(io::IO, ::$typ) = print(io, $body)
        parent(::$typ) = saturn
        naifid(::$typ) = $id
        export $sym, $typ
    end
end

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
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        Base.show(io::IO, ::$typ) = print(io, $body)
        parent(::$typ) = uranus
        naifid(::$typ) = $id
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
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        Base.show(io::IO, ::$typ) = print(io, $body)
        parent(::$typ) = neptune
        naifid(::$typ) = $id
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
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        Base.show(io::IO, ::$typ) = print(io, $body)
        parent(::$typ) = pluto
        naifid(::$typ) = $id
        export $sym, $typ
    end
end

