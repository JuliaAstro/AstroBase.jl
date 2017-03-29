struct RAC <: Frame end

function Rotation(::Type{GCRF}, ::Type{RAC}, r, v)
    normal = cross(r, v)
    normal = normal / norm(normal)
    tangential = v / norm(v)
    orthogonal = cross(v, normal)
    orthogonal = orthogonal / norm(orthogonal)
    m = hcat(orthogonal, tangential, normal)
    Rotation{GCRF, RAC}(m)
end
