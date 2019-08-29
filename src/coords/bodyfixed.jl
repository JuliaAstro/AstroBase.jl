export RAC

using ..Frames: AbstractFrame, ICRF

struct RAC <: AbstractFrame end

function Rotation(::RAC, ::ICRF, r, v)
    normal = cross(r, v)
    normal = normal / norm(normal)
    tangential = v / norm(v)
    orthogonal = cross(v, normal)
    orthogonal = orthogonal / norm(orthogonal)
    m = hcat(orthogonal, tangential, normal)
    Rotation{RAC(), icrf}(m)
end
