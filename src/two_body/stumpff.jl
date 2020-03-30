#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
using SpecialFunctions: gamma
export c2, c3

function c2(psi)
    eps = 1.0
    if psi > eps
        res = (1 - cos(sqrt(psi))) / psi
    elseif psi < -eps
        res = (cosh(sqrt(-psi)) - 1) / (-psi)
    else
        res = 1.0 / 2.0
        delta = (-psi) / gamma(2 + 2 + 1)
        k = 1
        while res + delta != res
            res += delta
            k += 1
            delta = (-psi)^k / gamma(2*k + 2 + 1)
        end
    end
    return res
end

function c3(psi)
    eps = 1.0
    if psi > eps
        res = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi))
    elseif psi < -eps
        res = (sinh(sqrt(-psi)) - sqrt(-psi)) / (-psi * sqrt(-psi))
    else
        res = 1.0 / 6.0
        delta = (-psi) / gamma(2 + 3 + 1)
        k = 1
        while res + delta != res
            res += delta
            k += 1
            delta = (-psi)^k / gamma(2*k + 3 + 1)
        end
    end
    return res
end
