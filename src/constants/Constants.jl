#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
module Constants

export astronomical_unit

const AU = 1.49597870700e8 # km

astronomical_unit(::Type{Float64}) = AU
astronomical_unit() = astronomical_unit(Float64)

end

