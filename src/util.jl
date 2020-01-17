module Util

using AccurateArithmetic: dot_oro
using LinearAlgebra: norm, ⋅, ×, normalize
using StaticArrays: SVector

export
    sec2rad, rad2sec, normalize_angle, angle,
    dms2rad, rad2dms, sec2deg, deg2sec,
    plane_section, point_on_limb

"""
    sec2rad(sec)

Convert an angle in arcseconds to radians.

# Example

```jldoctest
julia> sec2rad(3600 * 30)
0.5235987755982988
```
"""
sec2rad(sec) = deg2rad(sec / 3600)

"""
    rad2sec(rad)
Convert an angle in radians to arcseconds.

# Example

```jldoctest
julia> rad2sec(0.5235987755982988)
107999.99999999999
```
"""
rad2sec(rad) = rad2deg(rad) * 3600

sec2deg(sec) = sec / 3600
deg2sec(deg) = deg * 3600

function normalize_angle(angle, center)
    angle - 2π * floor((angle + π - center) / 2π)
end

function dms2rad(deg, arcmin, arcsec)
    deg2rad(deg + arcmin/60 + arcsec/3600)
end

function rad2dms(rad)
    d = rad2deg(rad)
    deg = trunc(d)
    arcmin = trunc((d-deg)*60)
    arcsec = (d-deg-arcmin/60)*3600
    return deg, arcmin, arcsec
end

function angle(v1, v2)
    normprod = norm(v1) * norm(v2)
    if normprod == 0.0
        throw(DomainError())
    else
        v1v2 = v1 ⋅ v2
        threshold = normprod * 0.9999
        if v1v2 >= -threshold && v1v2 <= threshold
            return acos(v1v2 / normprod)
        else
            v3n = norm(v1 × v2)
            return v1v2 >= 0.0 ? (v3n / normprod) : π - asin(v3n / normprod)
        end
    end
end

function azimuth(v)
    atan(v[2], v[1])
end

function elevation(v)
    asin(v[3] / norm(v))
end

function vector_azel(az, el)
    saz, caz = sincos(az)
    sel, cel = sincos(el)
    SVector(caz * cel, saz * cel, sel)
end

function angular_velocity(ψ, δψ, θ, δθ, ϕ, δϕ)
    Ω₁ = δψ * sin(θ) * sin(ϕ) - δθ * cos(ϕ)
    Ω₂ = δψ * sin(θ) * cos(ϕ) + δθ * sin(ϕ)
    Ω₃ = δψ * cos(θ) + δϕ
    [Ω₁, Ω₂, Ω₃]
end

@inbounds function frame(x)
    length(x) < 3 && throw(ArgumentError("Input vector must have three elements."))
    iszero(x) && return [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]

    x̂ = normalize(x[1:3])
    ŷ = Array{Float64}(undef, 3)
    ẑ = Array{Float64}(undef, 3)

    x2 = x̂ .^ 2
    _, idx = findmin(x2)
    s1, s2, s3 = circshift(1:3, 1-idx)
    f = sqrt(x2[s2] + x2[s3])

    ŷ[s1] = 0.0
    ŷ[s2] = -x̂[s3] / f
    ŷ[s3] = x̂[s2] / f
    ẑ[s1] = f
    ẑ[s2] = -x̂[s1] * ŷ[s3]
    ẑ[s3] = x̂[s1] * ŷ[s2]

    return x̂, ŷ, ẑ
end

function orthogonal(v)
    threshold = 0.6 * norm(v)
    iszero(threshold) && throw(ArgumentError("Norm of `v` is zero."))

    x, y, z = v
    if abs(x) <= threshold
        inverse = 1 / sqrt(y^2 + z^2)
        return [0.0, inverse * z, -inverse * y]
    elseif abs(y) <= threshold
        inverse = 1 / sqrt(x^2 + z^2)
        return [-inverse * z, 0.0, inverse * x]
    end
    inverse = 1 / sqrt(x^2 + y^2)
    return [inverse * y, -inverse * x, 0.0]
end

function plane_section(ellipsoid, plane_point, plane_normal)
    u = orthogonal(plane_normal)
    v = normalize(plane_normal × u)
    ue = u ./ ellipsoid
    ve = v ./ellipsoid
    pe = plane_point ./ ellipsoid
    alpha = ue ⋅ ue
    beta = ve ⋅ ve
    gamma = dot_oro(ue, ve)
    delta = dot_oro(pe, ue)
    epsilon = dot_oro(pe, ve)
    zeta = dot_oro([pe; 1.0], [pe; -1.0])

    if abs(gamma) < floatmin(Float64)
        tan_theta = 0.0
    else
        bma = beta - alpha
        tan_theta = bma >= 0.0 ? -2gamma / (bma + sqrt(bma^2 + 4gamma^2)) : -2gamma / (bma - sqrt(bma^2 + 4gamma^2))
    end
    tan2 = tan_theta^2
    cos2 = 1 / (1 + tan2)
    sin2 = tan2 * cos2
    cos_sin = tan_theta * cos2
    co = sqrt(cos2)
    si = tan_theta * co

    denom = dot_oro([gamma, -alpha], [gamma, beta])
    tau_c = dot_oro([beta, -gamma], [delta, epsilon]) / denom
    nu_c = dot_oro([alpha, -gamma], [epsilon, delta]) / denom

    twogcs = 2 * gamma * cos_sin
    A = alpha * cos2 + beta * sin2 + twogcs
    B = alpha * sin2 + beta * cos2 - twogcs
    F = (alpha * tau_c + 2 * (gamma * nu_c + delta)) * tau_c + (beta * nu_c + 2 * epsilon) * nu_c + zeta
    l = sqrt(-F / A)
    m = sqrt(-F / B)

    isnan(l + m) && return nothing

    if l > m
        return plane_point .+ tau_c .* u .+ nu_c .* v, co .* u .+ si .* v, -si .* u .+ co .* v, l, m
    else
        return plane_point .+ tau_c .* u .+ nu_c .* v, si .* u .- co .* v, co .* u .+ si .* v, m, l
    end
end

distance(v1, v2) = norm(v1 .- v2)

function isinside(ellipsoid, point)
    scaled = ellipsoid .- point
    return scaled ⋅ scaled <= 1.0
end

function point_on_limb(ellipsoid, observer, outside)
    isinside(ellipsoid, observer) && throw(ArgumentError("`observer` is inside `ellipsoid`."))

    normal = observer × outside
    center, u, v, a, b = plane_section(ellipsoid, zeros(3), normal)
    a2 = a^2
    b2 = b^2
    delta = observer - center
    xo = delta ⋅ u
    yo = delta ⋅ v
    xo2 = xo^2
    yo2 = yo^2
    alpha = a2 * yo2 + b2 * xo2
    beta = a2 * b2 * xo
    gamma = a2 * a2 * (b2 - yo2)
    sq = sqrt(beta^2 - alpha * gamma)
    if beta > 0
        s = beta + sq
        xt1 = s / alpha
        xt2 = gamma / s
    else
        s = beta - sq
        xt1 = gamma / s
        xt2 = s / alpha
    end
    t1 = center .+ xt1 .* u .+ (b2 * (a2 - xt1 * xo) / (a2 * yo)) .* v
    t2 = center .+ xt2 .* u .+ (b2 * (a2 - xt2 * xo) / (a2 * yo)) .* v
    return distance(t1, outside) <= distance(t2, outside) ? t1 : t2
end

end

