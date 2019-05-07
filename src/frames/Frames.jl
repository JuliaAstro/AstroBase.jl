module Frames

using ..Bodies:
    ALL_NAMES

using ItemGraphs: ItemGraph, add_edge!, items

export AbstractFrame, InertialFrame, RotatingFrame, ICRF, icrf, from_sym, @frame

abstract type AbstractFrame end
abstract type InertialFrame <: AbstractFrame end
abstract type RotatingFrame <: AbstractFrame end

Base.show(io::IO, frame::AbstractFrame) = print(io, string(nameof(typeof(frame))))

const FRAMES = ItemGraph{Symbol}()
from_sym(sym::Symbol) = from_sym(Val(sym))

struct ICRF <: InertialFrame end
const icrf = ICRF()
from_sym(::Val{:ICRF}) = icrf
path_frames(::F1, ::F2) where {F1, F2} = items(FRAMES, nameof(F1), nameof(F2))

macro frame(name::Symbol, cname::Symbol, parent::Symbol, typ::Symbol)
    quote
        struct $(esc(name)) <: $(esc(typ)) end
        const $(esc(cname)) = $(esc(name))()
        add_edge!(FRAMES, $(Meta.quot(name)), $(Meta.quot(parent)))
        add_edge!(FRAMES, $(Meta.quot(parent)), $(Meta.quot(name)))
        Frames.from_sym(::Val{$(Meta.quot(name))}) = $(esc(cname))
    end
end

include("rotations.jl")

for name in ALL_NAMES
    endswith(name, "Barycenter") && continue

    frame = Symbol("IAU", name)
    cname = Symbol("iau_", lowercase(name))
    @eval begin
        struct $frame <: RotatingFrame end
        const $cname = $frame()
        from_sym(::Val{$frame}) = $cname
        export $frame, $cname
    end
end

function __init__()
    for name in ALL_NAMES
        endswith(name, "Barycenter") && continue

        frame = Symbol("IAU", name)
        cname = Symbol("iau_", lowercase(name))
        add_edge!(FRAMES, frame, :ICRF)
        add_edge!(FRAMES, :ICRF, frame)
    end
end

end

