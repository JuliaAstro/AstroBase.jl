export keplerian, cartesian, perifocal, isprograde, isretrograde, ispolar

keplerian(rv, μ) = keplerian(rv[1:3], rv[4:6], μ)

function keplerian(r, v, µ)
    rm = norm(r)
    vm = norm(v)
    h = cross(r, v)
    hm = norm(h)
    k = [0.0, 0.0, 1.0]
    n = cross(k, h)
    nm = norm(n)
    xi = vm^2/2 - µ/rm
    ec = ((vm^2 - µ/rm)*r - v*dot(r, v))/µ
    ecc = norm(ec)
    inc = angle(h, k)

    equatorial = abs(inc) ≈ 0
    circular = ecc ≈ 0

    if circular
        # Semi-latus rectum
        sma = hm^2/µ
    else
        sma = -µ/(2*xi)
    end

    if equatorial && !circular
        node = 0.0
        # Longitude of pericenter
        peri = mod2pi(atan2(ec[2], ec[1]))
        ano = mod2pi(atan2(h⋅cross(ec, r) / hm, r⋅ec))
    elseif !equatorial && circular
        node = mod2pi(atan2(n[2], n[1]))
        peri = 0.0
        # Argument of latitude
        ano = mod2pi(atan2(r⋅cross(h, n) / hm, r⋅n))
    elseif equatorial && circular
        node = 0.0
        peri = 0.0
        # True longitude
        ano = mod2pi(atan2(r[2], r[1]))
    else
        node = mod2pi(atan2(n[2], n[1]))
        # peri = mod2pi(atan2(ec⋅cross(h, n) / hm, ec⋅n))
        # ano = mod2pi(atan2(r⋅cross(h, ec) / hm, r⋅ec))
        if sma > 0.0
            ese = (r ⋅ v) / sqrt(μ * sma)
            ece = rm * vm^2 / μ - 1.0
            ano = ecctotrue(atan2(ese, ece), ecc)
        else
            # TODO: Implement this.
            error("Hyperbolic elements not implemented.")
        end
        rnode = [cos(node), sin(node), 0.0]
        px = r ⋅ rnode
        py = r ⋅ (h × rnode) / hm
        peri = atan2(py, px) - ano
    end

    sma, ecc, inc, node, peri, ano
end

function cartesian(sma, ecc, inc, node, peri, ano, μ)
    if ecc ≈ 0
        p = sma
    else
        p = sma*(1 - ecc^2)
    end

    r_pqw, v_pqw = perifocal(p, ecc, ano, μ)
    M = rotation_matrix(313, -peri, -inc, -node)
    return M * r_pqw, M * v_pqw
end

function perifocal(p, ecc, ano, μ)
    r_pqw = [p * cos(ano) / (1 + ecc * cos(ano)), p * sin(ano) / (1 + ecc * cos(ano)), zero(p)]
    v_pqw = [-sqrt(μ/p) * sin(ano), sqrt(μ/p) * (ecc + cos(ano)), zero(sqrt(μ/p))]
    return r_pqw, v_pqw
end

cartesian(el, mu) = cartesian(el..., mu)

isprograde(inc) = inc < π/2 && !(inc ≈ π/2)
isretrograde(inc) = inc > π/2 && !(inc ≈ π/2)
ispolar(inc) = inc ≈ π/2
