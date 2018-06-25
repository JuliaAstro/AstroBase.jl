module AstroBase

using Rotations
export c2ixys
# package code goes here
function c2ixys(x, y, s)
    r2 = x^2 + y^2
    e = r2 > 0.0? atan2(y, x) : 0.0
    d = atan(sqrt(r2 / (1.0 - r2)))
    RotZYZ(e, d, -(e+s))
end

end # module
