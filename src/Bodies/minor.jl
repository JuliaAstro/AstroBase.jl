#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

@body pluto_barycenter 9 Barycenter parent=ssb _export=true
@body pluto 999 MinorBody parent=pluto_barycenter _export=true

# We have data for these minor bodies from pck00010.tpc and gm_de431.tpc.

const MINOR_BODY_NAMES = (
    "Borrelly",
    "Halley",
    "Tempel1",
    "Wild2",
    "Ceres",
    "Pallas",
    "Vesta",
    "Lutetia",
    "Kleopatra",
    "Mathilde",
    "Eros",
    "Davida",
    "Steins",
    "Toutatis",
    "Itokawa",
    "Ida",
    "Gaspra",
    "Psyche",
)

const MINOR_BODY_IDS = (
    1000005,
    1000036,
    1000093,
    1000107,
    2000001,
    2000002,
    2000004,
    2000021,
    2000216,
    2000253,
    2000433,
    2000511,
    2002867,
    2004179,
    2025143,
    2431010,
    9511010,
    2000016,
)

for (id, body) in zip(MINOR_BODY_IDS, MINOR_BODY_NAMES)
    name = Symbol(lowercase(body))
    @eval begin
        @body $name $id MinorBody parent=ssb _export=true
    end
end

