#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module Astrometry

using LinearAlgebra: ⋅, normalize
using ..Constants: schwarzschild_radius_sun

export aberration

"""
    aberration(pnat, v, s, bm1)

Apply aberration to transform natural direction into proper direction.

# Arguments

- `pnat`: Natural direction to the source (unit vector)
- `v`: Observer barycentric velocity in units of c
- `s`: Distance between the Sun and the observer (au)
- `bm1`: : Reciprocal of Lorenz factor ``\\sqrt{1-|v|^2}``

# Output

Returns the proper direction to source (unit vector).

# Example

```jldoctest
julia> pnat = [-0.76321968546737951,-0.60869453983060384,-0.21676408580639883];

julia> v = [2.1044018893653786e-5,-8.9108923304429319e-5,-3.8633714797716569e-5];

julia> s = 0.99980921395708788;

julia> bm1 = 0.99999999506209258;

julia> act = aberration(pnat, v, s, bm1)
3-element Array{Float64,1}:
 -0.7631631094219556
 -0.6087553082505591
 -0.21679262693684712
```

# References

- [ERFA - ab](https://github.com/liberfa/erfa/blob/master/src/ab.c)
"""
function aberration(pnat, v, s, bm1)
    pdv = pnat ⋅ v
    w1 = 1.0 + pdv / (1.0 + bm1)
    w2 = schwarzschild_radius_sun() / s
    p = @. pnat * bm1 + w1 * v + w2 * (v - pdv * pnat)
    return normalize(p)
end

end
