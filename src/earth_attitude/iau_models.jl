#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

export iau1980, iau1982, iau1994, iau2000, iau2000a, iau2000b, iau2006, iau2006a

abstract type IAUModel end
abstract type IAU1980Model <: IAUModel end
abstract type IAU2000Model <: IAUModel end
abstract type IAU2006Model <: IAUModel end

#############################
# IAU 1976/1980/1994 models #
#############################

struct IAU1980 <: IAU1980Model end

"""
    `iau1980`

The singleton instance of type `IAU1980`, representing the IAU 1980 family of models.
"""
const iau1980 = IAU1980()

struct IAU1982 <: IAU1980Model end

"""
    `iau1982`

The singleton instance of type `IAU1982`, representing the IAU 1982 family of models.
"""
const iau1982 = IAU1982()

struct IAU1994 <: IAU1980Model end

"""
    `iau1994`

The singleton instance of type `IAU1994`, representing the IAU 1994 family of models.
"""
const iau1994 = IAU1994()

###################
# IAU 2000 models #
###################

struct IAU2000 <: IAU2000Model end

"""
    `iau2000`

The singleton instance of type `IAU2000`, representing the IAU 2000 family of models.
"""
const iau2000 = IAU2000()

struct IAU2000A <: IAU2000Model end

"""
    `iau2000a`

The singleton instance of type `IAU2000A`, representing the IAU 2000A family of models.
"""
const iau2000a = IAU2000A()

struct IAU2000B <: IAU2000Model end

"""
    `iau2000b`

The singleton instance of type `IAU2000B`, representing the IAU 2000B family of models.
"""
const iau2000b = IAU2000B()

###################
# IAU 2006 models #
###################

struct IAU2006 <: IAU2006Model end

"""
    `iau2006`

The singleton instance of type `IAU2006`, representing the IAU 2006 family of models.
"""
const iau2006 = IAU2006()

struct IAU2006A <: IAU2006Model end

"""
    `iau2006a`

The singleton instance of type `IAU2006a`, representing the IAU 2006A family of models.
"""
const iau2006a = IAU2006A()

