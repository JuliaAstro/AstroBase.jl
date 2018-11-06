import Base: ∘, inv

export ComposedRotation, Rotation, ∘, origin, target

abstract type AbstractRotation end

struct Rotation{F1<:Frame,F2<:Frame} <: AbstractRotation
    m::Matrix
    δm::Matrix
end
Rotation{F1, F2}(m) where {F1<:Frame, F2<:Frame} = Rotation{F1, F2}(m, zeros(3,3))

function Rotation(::Type{F}, ::Type{F}, ep::Epoch) where F<:Frame
    Rotation{F,F}(eye(3), zeros(3,3))
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

function ∘(rot1::Rotation{F1,F}, rot2::Rotation{F,F2}) where {F<:Frame,F1<:Frame,F2<:Frame}
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

function haspath(graph, origin, target)
    haspath = false
    queue = [origin]
    links = Dict{DataType, DataType}()
    while !isempty(queue)
        node = shift!(queue)
        if node == target
            break
        end
        for neighbour in graph[node]
            if !haskey(links, neighbour)
                push!(queue, neighbour)
                merge!(links, Dict{DataType, DataType}(neighbour=>node))
            end
        end
    end
    if haskey(links, target)
        haspath = true
    end
    return haspath
end

function findpath(graph, origin, target)
    if isempty(graph[origin])
        error("There are no convert methods with source type '$origin' defined.")
    end
    if !haspath(graph, origin, target)
        error("No conversion path '$origin' -> '$target' found.")
    end
    queue = PriorityQueue(DataType, Int)
    prev = Dict{DataType,Nullable{DataType}}()
    distance = Dict{DataType, Int}()
    for node in keys(graph)
        merge!(prev, Dict(node=>Nullable{DataType}()))
        merge!(distance, Dict(node=>typemax(Int)))
        enqueue!(queue, node, distance[node])
    end
    distance[origin] = 0
    queue[origin] = 0
    while !isempty(queue)
        node = dequeue!(queue)
        node == target && break
        for neighbour in graph[node]
            alt = distance[node] + 1
            if alt < distance[neighbour]
                distance[neighbour] = alt
                prev[neighbour] = Nullable(node)
                queue[neighbour] = alt
            end
        end
    end
    path = DataType[]
    current = target
    while !isnull(prev[current])
        unshift!(path, current)
        current = get(prev[current])
    end
    return path
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

@generated function Rotation(::Type{F1}, ::Type{F2}, ep) where {F1<:Frame,F2<:Frame}
    gen_rotation(F1, F2, ep)
end
