export RAC

struct RAC <: Frame end

function Rotation(::Type{RAC}, ::Type{GCRF}, r, v)
    normal = cross(r, v)
    normal = normal / norm(normal)
    tangential = v / norm(v)
    orthogonal = cross(v, normal)
    orthogonal = orthogonal / norm(orthogonal)
    m = hcat(orthogonal, tangential, normal)
    Rotation{RAC, GCRF}(m)
end
