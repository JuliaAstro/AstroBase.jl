#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

const BODIES = ItemGraph{NAIFId, NAIFId}(SimpleGraph())

register_body!(id) = add_vertex!(BODIES, id)
link_bodies!(id1, id2) = add_edge!(BODIES, id1, id2)

path_ids(from::CelestialBody, to::CelestialBody) = items(BODIES, naifid(from), naifid(to))

"""
    @body(name, id, supertype, type=Name, parent=nothing)

Define a new celestial `body` with NAIF `id` and `supertype` which is the singleton
instance of `type`. Optionally, provide a `parent` body or pseudo-body.

# Example

```jldoctest
julia> @body rupert 42 MinorBody type=Persephone parent=ssb

julia> rupert isa Persephone
true

julia> typeof(rupert) <: MinorBody
true

julia> parent(rupert)
Solar System Barycenter

julia> naifid(rupert)
42
```
"""
macro body(name::Symbol, id::Int, super::Symbol, args...)
    name_comps = uppercasefirst.(split(String(name), '_'))
    str = join(name_comps, " ")
    typ = Symbol(join(name_comps))
    parent = nothing
    if !(super in (:CelestialBody, :Barycenter, :Planet, :NaturalSatellite, :MinorBody))
        throw(ArgumentError("Invalid supertype: $super"))
    end
    for a in args
        a isa Expr || continue
        a.head == :(=) || continue

        if a.args[1] == :type
            val = a.args[2]
            val isa Symbol || throw(ArgumentError("Invalid argument: $a"))
            typ = val
            str = join(split(String(typ), r"(?=[A-Z])"), " ")
        elseif a.args[1] == :parent
            val = a.args[2]
            val isa Symbol || throw(ArgumentError("Invalid argument: $a"))
            parent = val
        end
    end
    reg = if parent === nothing
        :(register_body!($id))
    else
        :(link_bodies!(naifid($parent), $id))
    end
    id_expr = :(naifid(::$typ) = $id)
    fromid_expr = :(from_naifid(::Val{$id}) = $name)
    return quote
        struct $(esc(typ)) <: $(esc(super)) end
        const $(esc(name)) = $(esc(typ))()
        Base.show(io::IO, ::$(esc(typ))) = print(io, "$($str)")
        $(esc(id_expr))
        $(esc(fromid_expr))
        $(esc(reg))
        nothing
    end
end

@body ssb 0 Barycenter type=SolarSystemBarycenter
@body sun 10 CelestialBody parent=ssb

# struct SolarSystemBarycenter <: Barycenter end
# const ssb = SolarSystemBarycenter()
# Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
# naifid(::SolarSystemBarycenter) = 0
# from_naifid(::Val{0}) = ssb
# add_vertex!(BODIES, 0)
#
# struct Sun <: CelestialBody end
# const sun = Sun()
# parent(::Sun) = ssb
# naifid(::Sun) = 10
# from_naifid(::Val{10}) = sun
# add_edge!(BODIES, 0, 10)

