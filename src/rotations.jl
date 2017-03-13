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
        m.sig isa UnionAll && continue
        length(m.sig.parameters) != 4 && continue

        t1 = m.sig.parameters[2].parameters[1]
        t2 = m.sig.parameters[3].parameters[1]
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

function gen_rotation(F1, F2, t)
    graph = getgraph()
    if !Convertible.haspath(graph, F1, F2)
        error("No conversion path '$F1' -> '$F2' found.")
    end
    path = Convertible.findpath(graph, F1, F2)
    ex = :(Rotation($F1, $(path[1]), t))
    for i in eachindex(path[2:end])
        t1 = path[i]
        t2 = path[i+1]
        ex = :(compose($ex, Rotation($t1, $t2, t)))
    end
    return ex
end

@generated function Rotation{F1<:Frame,F2<:Frame}(::Type{F1}, ::Type{F2}, t::Float64)
    gen_transform(F1, F2, t)
end
