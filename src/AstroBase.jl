module AstroBase

using Rotations
export pom00
@inline function pom00(rx, ry, sp)
    RotZYX(sp, -rx, -ry)
end

end
