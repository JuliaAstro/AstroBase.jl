import Convertible: haspath, findpath
import Base: ∘, inv

export ComposedRotation, Rotation, ∘, origin, target

abstract type AbstractRotation end

struct Rotation{F1<:Frame,F2<:Frame} <: AbstractRotation
    m::Matrix{Float64}
    δm::Matrix{typeof(1.0/s)}
end
Rotation{F1, F2}(m) where {F1<:Frame, F2<:Frame} = Rotation{F1, F2}(m, zeros(3,3)*(1/s))

function Rotation(::Type{F}, ::Type{F}, ep::Epoch) where F<:Frame
    Rotation{F,F}(eye(3), zeros(3,3)*(1/s))
end

origin(::Rotation{F1,<:Frame}) where F1<:Frame = F1
target(::Rotation{<:Frame,F2}) where F2<:Frame = F2

function inv(rot::Rotation{F1, F2}) where {F1<:Frame, F2<:Frame}
    Rotation{F2, F1}(rot.m', rot.δm')
end

(rot::Rotation)(r, v) = rot.m * r, rot.δm * r + rot.m * v
(rot::Rotation)(rv) = rot.m * rv[1], rot.δm * rv[1] + rot.m * rv[2]

struct ComposedRotation{R1<:AbstractRotation,R2<:AbstractRotation} <: AbstractRotation
    rot1::R1
    rot2::R2
end
origin(rot::ComposedRotation) = origin(rot.rot1)
target(rot::ComposedRotation) = target(rot.rot2)

(rot::ComposedRotation)(args...) = rot.rot2(rot.rot1(args...))

inv(rot::ComposedRotation) = ComposedRotation(inv(rot.rot2), inv(rot.rot1))

function ∘{F<:Frame,F1<:Frame,F2<:Frame}(rot1::Rotation{F1,F}, rot2::Rotation{F,F2})
    ComposedRotation(rot1, rot2)
end

function getgraph()
    graph = Dict{DataType,Set{DataType}}()
    for m in methods(Rotation)
        params = m.sig isa UnionAll ? Base.unwrap_unionall(m.sig).parameters :
            m.sig.parameters
        length(params) != 4 && continue

        t1 = params[2].parameters[1]
        t2 = params[3].parameters[1]
        t1 = isleaftype(t1) ? t1 : t1.ub
        t2 = isleaftype(t2) ? t2 : t2.ub

        if !haskey(graph, t1)
            merge!(graph, Dict(t1=>Set{DataType}()))
        end
        if !haskey(graph, t2)
            merge!(graph, Dict(t2=>Set{DataType}()))
        end
        push!(graph[t1], t2)
    end
    graph
end

function gen_rotation(F1, F2, ep)
    graph = getgraph()
    FF1 = supertype(F1) == Frame ? F1 : supertype(F1)
    FF2 = supertype(F2) == Frame ? F2 : supertype(F2)
    if !haspath(graph, FF1, FF2)
        error("No conversion path '$F1' -> '$F2' found.")
    end
    path = findpath(graph, FF1, FF2)
    ex = :(Rotation($F1, $(path[1]), ep))
    for i in eachindex(path[2:end-1])
        t1 = path[i]
        t2 = path[i+1]
        ex = :(ComposedRotation($ex, Rotation($t1, $t2, ep)))
    end
    ex = :(ComposedRotation($ex, Rotation($(path[end-1]), $F2, ep)))
    return ex
end

@generated function Rotation{F1<:Frame,F2<:Frame}(::Type{F1}, ::Type{F2}, ep)
    gen_rotation(F1, F2, ep)
end
