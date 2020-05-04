#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

import Base: position

export AbstractEphemeris,
    position,
    position!,
    velocity,
    velocity!,
    state,
    state!

"""
    AbstractEphemeris

Abstract supertype for ephemerides.
"""
abstract type AbstractEphemeris end

function position! end
function velocity end
function velocity! end
function state end
function state! end

