#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
struct CubicSpline{XT,T}
    n::Int
    x::Vector{XT}
    y::Vector{T}
    c::Matrix{T}
end

"""
    CubicSpline(x, y)

Construct a cubic spline interpolator. `x` and `y` must have the same length and
at least four elements.

# Example

```jldoctest
julia> x = 1:5;

julia> y = [2.0, 5.0, 3.2, 4.8, 6.7];

julia> spl = CubicSpline(x, y);

julia> spl(3.5)
3.5921875
```
"""
function CubicSpline(x, y)
    n = length(x)
    n < 4 && throw(ArgumentError("At least four points are needed."))
    length(y) != n && throw(ArgumentError("`x` and `y` must have the same length."))

    dl = similar(y, n-1)
    d = similar(y, n)
    du = similar(y, n-1)
    b = similar(y, n)
    c = similar(y, 4, n-1)
    @views begin
        dx = diff(x)
        slope = diff(y) ./ dx

        d[2:end-1] .= 2 .* (dx[1:end-1] .+ dx[2:end])
        du[2:end] .= dx[1:end-1]
        dl[1:end-1] .= dx[2:end]
        b[2:end-1] .= 3 * (dx[2:end] .* slope[1:end-1]
            .+ dx[1:end-1] .* slope[2:end])

        # Not-a-knot boundary conditions
        d[1] = dx[2]
        du[1] = x[3] - x[1]
        δ = x[3] - x[1]
        b[1] = ((dx[1] + 2δ) * dx[2] * slope[1]
            + dx[1]^2 * slope[2]) / δ
        d[end] = dx[end-1]
        dl[end] = x[end] - x[end-2]
        δ = x[end] - x[end-2]
        b[end] = ((dx[end]^2 * slope[end-1]
                + (2δ + dx[end]) * dx[end-1] * slope[end]) / δ)

        A = Tridiagonal(dl, d, du)
        s = A \ b
        t = (s[1:end-1] .+ s[2:end] .- 2 .* slope) ./ dx

        c[1,:] .= y[1:end-1]
        c[2,:] .= s[1:end-1]
        c[3,:] .= (slope .- s[1:end-1]) ./ dx .- t
        c[4,:] .= t ./ dx
    end

    return CubicSpline(n, collect(x), y, c)
end

@inbounds function evaluate(spl::CubicSpline, x1)
    ascending = spl.x[end] >= spl.x[1]
    idx = if ascending
        searchsortedlast(spl.x, x1)
    else
        searchsortedlast(spl.x, x1, lt=!isless)
    end
    idx == 0 && return spl.y[1]
    idx == spl.n && return spl.y[end]
    return @evalpoly(x1 - spl.x[idx],
        spl.c[1,idx], spl.c[2,idx], spl.c[3,idx], spl.c[4,idx])
end

(spl::CubicSpline)(x1) = evaluate(spl, x1)

