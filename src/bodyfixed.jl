struct RAC <: Frame end

function Rotation(::Type{GCRF}, ::Type{RAC}, r, v)
    normal = cross(r, v)
    normal = normal / norm(normal)
    tangential = v / norm(v)
    orthogonal = cross(v, normal)
    orthogonal = orthogonal / norm(orthogonal)
    m = hcat(orthogonal, tangential, normal)
    δm = zeros(3, 3) * (1/s)
    Rotation{GCRF, RAC}(m, δm)
end
