abstract type Planet <: CelestialBody end

const PLANET_NAMES = (
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
)

for (i, body) in enumerate(PLANET_NAMES)
    typ = Symbol(body)
    sym = Symbol(lowercase(body))
    id = 100i + 99
    typ_bc = Symbol(body, "Barycenter")
    nam_bc = string(Symbol(body, " Barycenter"))
    sym_bc = Symbol(lowercase(body), "_barycenter")
    @eval begin
        struct $typ_bc <: Barycenter end
        const $sym_bc = $typ_bc()
        Base.show(io::IO, ::$typ_bc) = print(io, $nam_bc)
        parent(::$typ_bc) = ssb
        naifid(::$typ_bc) = $i
        from_naifid(::Val{$i}) = $sym_bc
        export $sym_bc, $typ_bc

        struct $typ <: Planet end
        const $sym = $typ()
        parent(::$typ) = $sym_bc
        naifid(::$typ) = $id
        from_naifid(::Val{$id}) = $sym
        export $sym, $typ
    end
end

# FIXME: Move this to a sensible place
export j2, mean_radius, polar_radius, equatorial_radius
j2(::Earth) = 1.08262668e-3
mean_radius(::Earth) = 6371.0084
equatorial_radius(::Earth) = 6378.1366
polar_radius(::Earth) = 6356.7519
