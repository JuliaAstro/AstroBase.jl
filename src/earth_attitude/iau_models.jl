export iau1980, iau2000a, iau2000b, iau2006

abstract type IAUModel end

struct IAU1980 <: IAUModel end
const iau1980 = IAU1980()

struct IAU2000A <: IAUModel end
const iau2000a = IAU2000A()

struct IAU2000B <: IAUModel end
const iau2000b = IAU2000B()

struct IAU2006 <: IAUModel end
const iau2006 = IAU2006()
