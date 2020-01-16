module Util

using LinearAlgebra: norm, ⋅, ×, normalize
using StaticArrays: SVector

export
    sec2rad, rad2sec, normalize_angle, angle,
    dms2rad, rad2dms, sec2deg, deg2sec

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

end

