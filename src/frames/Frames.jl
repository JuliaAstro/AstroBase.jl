module Frames

using ..Bodies:
    ALL_NAMES

using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, items

export AbstractFrame, InertialFrame, RotatingFrame, ICRF, icrf, from_sym, @frame

abstract type AbstractFrame end
abstract type InertialFrame <: AbstractFrame end
abstract type RotatingFrame <: AbstractFrame end

Base.show(io::IO, frame::AbstractFrame) = print(io, string(nameof(typeof(frame))))

const FRAMES = ItemGraph{Symbol}(SimpleGraph())
from_sym(sym::Symbol) = from_sym(Val(sym))

struct ICRF <: InertialFrame end
const icrf = ICRF()
from_sym(::Val{:ICRF}) = icrf

function path_frames(frame1, frame2)
    name1 = nameof(typeof(frame1))
    name2 = nameof(typeof(frame2))
    items(FRAMES, name1, name2)
end

macro frame(name::Symbol, cname::Symbol, parent::Symbol, typ::Symbol)
    quote
        struct $(esc(name)) <: $(esc(typ)) end
        const $(esc(cname)) = $(esc(name))()
        add_edge!(FRAMES, $(Meta.quot(name)), $(Meta.quot(parent)))
        Frames.from_sym(::Val{$(Meta.quot(name))}) = $(esc(cname))
    end
end

include("rotations.jl")
include("iau.jl")
include("iers.jl")

end

