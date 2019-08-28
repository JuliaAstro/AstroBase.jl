using ..Frames: Rotation

# F1 -> F2
function transform(ep, rv, _, ::S, ::F1, ::C, ::S, ::F2, ::C) where {S, F1, F2, C}
    rot = Rotation(F1(), F2(), ep)
    State(ep, rot(rv)..., frame=F2(), body=C())
end

# S1 -> S2
function transform(ep, rv, _, ::S1, ::F, ::C, ::S2, ::F, ::C) where {S1, S2, F, C}
    State(Epoch{S2()}(ep), rv..., frame=F(), body=C())
end

# C1 -> C2
function transform(ep, rv, eph, ::S, ::F, ::C1, ::S, ::F, ::C2) where {S, F, C1, C2}
    rot = Rotation(F(), icrf, ep)
    rv_icrf = rot(rv)
    rv_body = state(eph, ep, C1(), C2())
    rv′ = inv(rot)(rv_icrf .- rv_body)
    State(ep, rv′..., frame=F(), body=C2())
end

# F1 -> F2, S1 -> S2
function transform(ep, rv, _, ::S1, ::F1, ::C, ::S2, ::F2, ::C) where {S1, S2, F1, F2, C}
    rot = Rotation(F1(), F2(), ep)
    State(Epoch{S2()}(ep), rot(rv)..., frame=F2(), body=C())
end

# F1 -> F2, C1 -> C2
function transform(ep, rv, eph, ::S, ::F1, ::C1, ::S, ::F2, ::C2) where {S, F1, F2, C1, C2}
    rot1 = Rotation(F1(), icrf, ep)
    rot2 = Rotation(icrf, F2(), ep)
    rv_icrf = rot1(rv)
    rv_body = state(eph, ep, C1(), C2())
    rv′ = rot2(rv_icrf .- rv_body)
    State(ep, rv′..., frame=F2(), body=C2())
end

# S1 -> S2, C1 -> C2
function transform(ep, rv, eph, ::S1, ::F, ::C1, ::S2, ::F, ::C2) where {S1, S2, F, C1, C2}
    rot = Rotation(F(), icrf, ep)
    rv_icrf = rot(rv)
    rv_body = state(eph, ep, C1(), C2())
    rv′ = inv(rot)(rv_icrf .- rv_body)
    State(Epoch{S2()}(ep), rv′..., frame=F(), body=C2())
end

# F1 -> F2, S1 -> S2, C1 -> C2
function transform(ep, rv, eph, ::S1, ::F1, ::C1, ::S2, ::F2, ::C2) where {S1, S2, F1, F2, C1, C2}
    rot1 = Rotation(F1(), icrf, ep)
    rot2 = Rotation(icrf, F2(), ep)
    rv_icrf = rot1(rv)
    rv_body = state(eph, ep, C1(), C2())
    rv′ = rot2(rv_icrf .- rv_body)
    State(Epoch{S2()}(ep), rv′..., frame=F2(), body=C2())
end

