import Convertible: haspath, findpath

export Rotation

struct Rotation{F1<:Frame,F2<:Frame}
    matrix::Matrix{Float64}
end

(t::Rotation)(rv) = t.matrix*rv

function compose{F<:Frame,F1<:Frame,F2<:Frame}(t1::Rotation{F1,F}, t2::Rotation{F,F2})
    Rotation{F1,F2}(t2.matrix*t1.matrix)
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
        ex = :(compose($ex, Rotation($t1, $t2, ep)))
    end
    ex = :(compose($ex, Rotation($(path[end-1]), $F2, ep)))
    return ex
end

@generated function Rotation{F1<:Frame,F2<:Frame}(::Type{F1}, ::Type{F2}, ep)
    gen_rotation(F1, F2, ep)
end
