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
    name_str = String(name)
    typ_str = join(uppercasefirst.(split(name_str, '_')))
    typ = Symbol(typ_str)
    parent = nothing
    _export = false
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
        elseif a.args[1] == :parent
            val = a.args[2]
            val isa Symbol || throw(ArgumentError("Invalid argument: $a"))
            parent = val
        elseif a.args[1] == :_export
            val = a.args[2]
            val isa Bool || throw(ArgumentError("Invalid argument: $a"))
            _export = val
        end
    end
    parts = split(String(typ), r"(?=[A-Z1-9])")
    show_str = join(parts, " ")
    doc_str = parts[end] == "Barycenter" ? join(["the"; parts], " ") : show_str
    super_str = String(super)
    reg = if parent === nothing
        :(register_body!($id))
    else
        :(link_bodies!(naifid($parent), $id))
    end
    exp = _export ? :(export $name, $typ) : :()
    id_expr = :(naifid(::$typ) = $id)
    fromid_expr = :(from_naifid(::Val{$id}) = $name)
    return quote
        """
            $($typ_str) <: $($super_str)

        A type representing $($doc_str).
        """
        struct $(esc(typ)) <: $(esc(super)) end

        """
            $($name_str)

        The singleton instance of the [`$($typ_str)`](@ref) type.
        """
        const $(esc(name)) = $(esc(typ))()
        Base.show(io::IO, ::$(esc(typ))) = print(io, "$($show_str)")
        $(esc(id_expr))
        $(esc(fromid_expr))
        $(esc(reg))
        $(esc(exp))
        nothing
    end
end

@body ssb 0 Barycenter type=SolarSystemBarycenter _export=true
@body sun 10 CelestialBody parent=ssb _export=true

