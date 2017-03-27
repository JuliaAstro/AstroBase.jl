function period(s::State)
    ele = keplerian(s)
    period(ele[1], μ(body(s)))
end

keplerian(s::State) = keplerian(rv_array(s), μ(body(s)))
semimajoraxis(s::State) = keplerian(s)[1]
eccentricity(s::State) = keplerian(s)[2]
inclination(s::State) = keplerian(s)[3]
ascendingnode(s::State) = keplerian(s)[4]
argumentofpericenter(s::State) = keplerian(s)[5]
trueanomaly(s::State) = keplerian(s)[6]

rotation_matrix{F<:Frame}(::Type{F}, ::Type{F}, ep::Epoch) = eye(6, 6)

rotation_matrix(f1::Type{RAC}, f2::Type{GCRF}, y::Vector{Float64}) = rotation_matrix(f1, f2, y[1:3], y[4:6])
function rotation_matrix(::Type{RAC}, ::Type{GCRF}, r::Vector{Float64}, v::Vector{Float64})
    normal = cross(r, v)
    normal = normal / norm(normal)
    tangential = v / norm(v)
    orthogonal = cross(v, normal)
    orthogonal = orthogonal / norm(orthogonal)
    hcat(orthogonal, tangential, normal)
end
