#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

@body pluto_barycenter 9 Barycenter parent=ssb _export=true
@body pluto 999 MinorBody parent=pluto_barycenter _export=true

const MINOR_BODY_NAMES = (
    "Ceres",
    "Pallas",
    "Vesta",
    "Psyche",
    "Lutetia",
    "Ida",
    "Eros",
    "Davida",
    "Gaspra",
    "Steins",
    "Itokawa",
    "Tempel1",
    "Borrelly",
)

const MINOR_BODY_IDS = (
    2000001,
    2000002,
    2000004,
    2000016,
    2000021,
    2431010,
    2000433,
    2000511,
    9511010,
    2002867,
    2025143,
    1000093,
    1000005,
)

for (id, body) in zip(MINOR_BODY_IDS, MINOR_BODY_NAMES)
    name = Symbol(lowercase(body))
    @eval begin
        @body $name $id MinorBody parent=ssb _export=true
    end
end

