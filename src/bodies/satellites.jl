abstract type NaturalSatellite <: CelestialBody end

export deimos, Deimos, phobos, Phobos, luna, Luna, moon, charon, Charon

# Earth

struct Luna <: NaturalSatellite end
const luna = Luna()
const moon = luna
parent(::Luna) = earth_barycenter
naifid(::Luna) = 301
from_naifid(::Val{301}) = luna

# Mars

struct Phobos <: NaturalSatellite end
const phobos = Phobos()
parent(::Phobos) = mars_barycenter
naifid(::Phobos) = 401
from_naifid(::Val{401}) = phobos

struct Deimos <: NaturalSatellite end
const deimos = Deimos()
parent(::Deimos) = mars_barycenter
naifid(::Deimos) = 402
from_naifid(::Val{402}) = deimos

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
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = saturn_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
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
    @eval begin
        struct $typ <: NaturalSatellite end
        const $sym = $typ()
        parent(::$typ) = pluto_barycenter
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end
