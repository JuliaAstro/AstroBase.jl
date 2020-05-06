#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

const PLANET_NAMES = (
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
)

for (i, body) in enumerate(PLANET_NAMES)
    planet = Symbol(lowercase(body))
    barycenter = Symbol(planet, "_barycenter")
    id = 100i + 99
    @eval begin
        @body $barycenter $i Barycenter parent=ssb _export=true
        @body $planet $id Planet parent=$barycenter _export=true
    end
end

