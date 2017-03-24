immutable State{F<:Frame, T<:Timescale, C<:CelestialBody} <: AbstractState
    epoch::Epoch{T}
    rv::Vector{Float64}
    frame::Type{F}
    body::Type{C}
end

function show(io::IO, s::State)
    println(io, "State{$(s.frame),$(s.epoch.scale),$(s.body)}:")
    println(io, " Epoch: $(s.epoch)")
    println(io, " Frame: $(s.frame)")
    println(io, " Body:  $(s.body)\n")
    println(io, " x: $(s.rv[1])")
    println(io, " y: $(s.rv[2])")
    println(io, " z: $(s.rv[3])")
    println(io, " u: $(s.rv[4])")
    println(io, " v: $(s.rv[5])")
    println(io, " w: $(s.rv[6])\n")
    sma, ecc, inc, node, peri, ano = keplerian(s)
    println(io, " a: $sma")
    println(io, " e: $ecc")
    println(io, " i: $(rad2deg(inc))")
    println(io, " ω: $(rad2deg(node))")
    println(io, " Ω: $(rad2deg(peri))")
    print(io, " ν: $(rad2deg(ano))")
end

function State{F<:Frame, T<:Timescale, C<:CelestialBody}(ep::Epoch{T}, rv, frame::Type{F}=GCRF, body::Type{C}=Earth)
    State(ep, rv, frame, body)
end

function State{F<:Frame, T<:Timescale, C<:CelestialBody}(ep::Epoch{T}, r, v, frame::Type{F}=GCRF, body::Type{C}=Earth)
    State(ep, [r; v], frame, body)
end

function State{F<:Frame, T<:Timescale, C<:CelestialBody}(ep::Epoch{T}, x, y, z, vx, vy, vz, frame::Type{F}=GCRF, body::Type{C}=Earth)
    State(ep, [x, y, z, vx, vy, vz], frame, body)
end

function (==){F<:Frame, T<:Timescale, C<:CelestialBody}(a::State{F,T,C}, b::State{F,T,C})
    return a.epoch == b.epoch && a.rv == b.rv && a.frame == b.frame && a.body == a.body
end

function isapprox{F<:Frame, T<:Timescale, C<:CelestialBody}(a::State{F,T,C}, b::State{F,T,C})
    return a.epoch ≈ b.epoch && a.rv ≈ b.rv && a.frame == b.frame && a.body == a.body
end

body(s::State) = constants(s.body)
μ(s::State) = μ(constants(s.body))
rv_array(s::State) = s.rv
epoch(s::State) = s.epoch
reference_frame(s::State) = s.frame

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

function keplerian2state(ep, sma, ecc, inc, node, peri, ano, frame=GCRF, body=Earth)
    State(ep, cartesian(sma, ecc, inc, node, peri, ano, μ(body))..., frame, body)
end

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
