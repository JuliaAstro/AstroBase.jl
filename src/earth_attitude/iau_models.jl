export iau1980, iau2000a, iau2000b, iau2006

abstract type IAUModel end

struct IAU1980 <: IAUModel end

"""
    `iau1980`

The singleton instance of type `IAU1980`, representing the IAU 1980 family of models.
"""
const iau1980 = IAU1980()

struct IAU2000A <: IAUModel end

"""
    `iau2000a`

The singleton instance of type `IAU2000A`, representing the IAU 2000A family of models.
"""
const iau2000a = IAU2000A()

struct IAU2000B <: IAUModel end

"""
    `iau2000b`

The singleton instance of type `IAU2000B`, representing the IAU 2000B family of models.
"""
const iau2000b = IAU2000B()

struct IAU2006 <: IAUModel end

"""
    `iau2006`

The singleton instance of type `IAU2006`, representing the IAU 2006 family of models.
"""
const iau2006 = IAU2006()
