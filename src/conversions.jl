import Base: convert

convert(::Type{State{F, T, C}}, s::State{F, T, C}) where {
    F<:Frame, T<:Timescale, C<:CelestialBody} = s

# F1 -> F2
function convert(::Type{State{F2, T, C}}, s::State{F1, T, C}) where {
    F1<:Frame, F2<:Frame, T<:Timescale, C<:CelestialBody}
    rot = Rotation(F1, F2, epoch(s))
    rv1 = rot(rv(s))
    State(epoch(s), rv1[1:3], rv1[4:6], F2, body(s))
end

# T1 -> T2
function convert(::Type{State{F, T2, C}}, s::State{F, T1, C}) where {
    F<:Frame, T1<:Timescale, T2<:Timescale, C<:CelestialBody}
    State(Epoch(T2, s.epoch), s.rv, F, s.body)
end

# C1 -> C2
function convert(::Type{State{F, T, C2}}, s::State{F, T, C1}) where {
    F<:Frame, T<:Timescale, C1<:CelestialBody, C2<:CelestialBody}
    M1 = rotation_matrix(GCRF, F, s.epoch)
    M2 = rotation_matrix(F, GCRF, s.epoch)
    rv = M1*s.rv
    body1 = state(C1, s.epoch)
    body2 = state(C2, s.epoch)
    State(s.epoch, M2*(rv + body1 - body2), F, C2)
end

# F1 -> F2, T1 -> T2
function convert(::Type{State{F2, T2, C}}, s::State{F1, T1, C}) where {
    F1<:Frame, F2<:Frame, T1<:Timescale, T2<:Timescale, C<:CelestialBody}
    M = rotation_matrix(F2, F1, s.epoch)
    State(Epoch(T2, s.epoch), M * s.rv, F2, s.body)
end

# F1 -> F2, C1 -> C2
function convert(::Type{State{F2, T, C2}}, s::State{F1, T, C1}) where {
    F1<:Frame, F2<:Frame, T<:Timescale, C1<:CelestialBody, C2<:CelestialBody}
    M1 = rotation_matrix(GCRF, F1, s.epoch)
    M2 = rotation_matrix(F2, GCRF, s.epoch)
    rv = M1*s.rv
    body1 = state(C1, s.epoch)
    body2 = state(C2, s.epoch)
    State(s.epoch, M2*(rv + body1 - body2), F2, C2)
end

# T1 -> T2, C1 -> C2
function convert(::Type{State{F, T2, C2}}, s::State{F, T1, C1}) where {
    F<:Frame, T1<:Timescale, T2<:Timescale, C1<:CelestialBody, C2<:CelestialBody}
    M = rotation_matrix(F, GCRF, s.epoch)
    body1 = state(C1, s.epoch)
    body2 = state(C2, s.epoch)
    State(Epoch(T2, s.epoch), s.rv + M*body1 - M*body2, F, C2)
end

# F1 -> F2, T1 -> T2, C1 -> C2
function convert(::Type{State{F2,T2,C2}}, s::State{F1,T1,C1}) where {
    F1<:Frame, F2<:Frame, T1<:Timescale, T2<:Timescale,
    C1<:CelestialBody, C2<:CelestialBody}
    M1 = rotation_matrix(GCRF, F1, s.epoch)
    M2 = rotation_matrix(F2, GCRF, s.epoch)
    rv = M1*s.rv
    body1 = state(C1, s.epoch)
    body2 = state(C2, s.epoch)
    State(Epoch(T2, s.epoch), M2*(rv + body1 - body2), F2, C2)
end
