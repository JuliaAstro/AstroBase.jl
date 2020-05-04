#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, add_vertex!, items

export ICRF, icrf, @frame

const FRAMES = ItemGraph{AbstractFrame}(SimpleGraph())

path_frames(from, to) = items(FRAMES, from, to)

register_frame!(frame) = add_vertex!(FRAMES, frame)
link_frames!(frame1, frame2) = add_edge!(FRAMES, frame1, frame2)

"""
    @frame(name, type=NAMEFrame, parent=nothing, rotating=false)

Define a new reference frame `name` which is the singleton instance of `type`.
Optionally provide a `parent` frame and indicate whether the frame is `rotating`.

# Example

```jldoctest
julia> @frame inertial

julia> isinertial(inertial)
true

julia> typeof(inertial)
INERTIALFrame

julia> @frame rotating type=Rotating parent=inertial rotating=true

julia> isrotating(rotating)
true

julia> typeof(rotating)
Rotating
```
"""
macro frame(name, args...)
    typ = Symbol(uppercase(String(name)), "Frame")
    parent = nothing
    rotating = false
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
        elseif a.args[1] == :rotating
            val = a.args[2]
            val isa Bool || throw(ArgumentError("Invalid argument: $a"))
            rotating = val
        end
    end
    atype = rotating ? :RotatingFrame : :InertialFrame
    vexpr = :(register_frame!($(esc(name))))
    eexpr = :(link_frames!($(esc(parent)), $(esc(name))))
    reg = parent === nothing ? vexpr : eexpr
    return quote
        struct $(esc(typ)) <: $(esc(atype)) end
        const $(esc(name)) = $(esc(typ))()
        $reg
        nothing
    end
end

@frame icrf type=ICRF

"""
    ICRF

A type representing the International Celestial Reference Frame (ICRF).

# Reference

- [SOFA](http://www.iausofa.org/publications/sofa_pn.pdf)
"""
ICRF

"""
    icrf

The singleton instance of the [`ICRF`](@ref) type representing the International Celestial
Reference Frame (ICRF).

# Reference

- [SOFA](http://www.iausofa.org/publications/sofa_pn.pdf)
"""
icrf

