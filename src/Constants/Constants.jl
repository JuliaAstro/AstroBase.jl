#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module Constants

using ..Time:
    TCB, TCG, TDB, TT,
    BarycentricCoordinateTime,
    BarycentricDynamicalTime,
    GeocentricCoordinateTime,
    TerrestrialTime

export
    KiloMeter,
    Meter,
    astronomical_unit,
    drift_rate,
    earth_rotation_angle_j2000,
    earth_rotation_rate,
    equatorial_radius_earth,
    gaussian_gravitational_const,
    geocentric_gravitational_const,
    geoid_potential,
    gravitational_const,
    heliocentric_gravitational_const,
    j2_factor_earth,
    j2_rate_earth,
    mass_ratio_moon,
    mean_angular_velocity_earth,
    obliquity_j2000,
    speed_of_light,
    tdb0

"""
    Meter

A singleton type for returning constants in SI units.
"""
struct Meter end

"""
    KiloMeter

A singleton type for returning constants in km-based units (the default).
"""
struct KiloMeter end

const SPEED_OF_LIGHT = 2.99792458e8 # m/s

"""
    speed_of_light([[T=Float64], unit=KiloMeter()])

Return the speed of light ``c``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
speed_of_light(::Type{Float64}, ::Meter) = SPEED_OF_LIGHT
speed_of_light(::Type{Float64}, ::KiloMeter) = speed_of_light(Float64, Meter()) * 1e-3
speed_of_light() = speed_of_light(Float64, KiloMeter())

const GAUSSIAN_GRAVITATIONAL_CONSTANT = 1.720209895e-2 # dimensionless

"""
    gaussian_gravitational_const([T=Float64])

Return the Gaussian gravitational constant ``k``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
gaussian_gravitational_const(::Type{Float64}) = GAUSSIAN_GRAVITATIONAL_CONSTANT
gaussian_gravitational_const() = gaussian_gravitational_const(Float64)

const LG = 6.969290134e-10 # dimensionless

"""
    drift_rate([T=Float64,] TT, TCG)

Return the rate of drift between TT and TCG ``L_G``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
drift_rate(::Type{Float64}, ::TerrestrialTime, ::GeocentricCoordinateTime) = LG
drift_rate(::TerrestrialTime, ::GeocentricCoordinateTime) = drift_rate(Float64, TT, TCG)

const LB = 1.550519768e-8 # dimensionless

"""
    drift_rate([T=Float64,] TDB, TCB)

Return the rate of drift between TDB and TCB ``L_B``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
drift_rate(
    ::Type{Float64},
    ::BarycentricDynamicalTime,
    ::BarycentricCoordinateTime,
) = LB
drift_rate(
    ::BarycentricDynamicalTime,
    ::BarycentricCoordinateTime,
) = drift_rate(Float64, TDB, TCB)

const TDB0 = -6.55e-5 # s

"""
    tdb0([T=Float64])

Return TDB–TCB at T0.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
tdb0(::Type{Float64}) = TDB0
tdb0() = tdb0(Float64)

const THETA0 = 0.7790572732640 # revolutions

"""
    earth_rotation_angle_j2000([T=Float64])

Return the Earth rotation angle at J2000.0 ``\\theta_0``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
earth_rotation_angle_j2000(::Type{Float64}) = THETA0
earth_rotation_angle_j2000() = earth_rotation_angle_j2000(Float64)

const DELTA_THETA = 1.00273781191135448 # revolutions/UT1 day

"""
    earth_rotation_rate([T=Float64])

Return the rate of advance of Earth rotation angle ``d\\theta/dUT1``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
earth_rotation_rate(::Type{Float64}) = DELTA_THETA
earth_rotation_rate() = earth_rotation_rate(Float64)

const G = 6.67428e-11 # m³/(kg * s²)
const G_ERR = 6.7e-15 # m³/(kg * s²)

"""
    gravitational_const([[T=Float64], unit=KiloMeter()])

Return the gravitational constant ``G`` and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
gravitational_const(::Type{Float64}, ::Meter) = G, G_ERR
gravitational_const(::Type{Float64}, ::KiloMeter) = gravitational_const(Float64, Meter()) .* 1e-9
gravitational_const() = gravitational_const(Float64, KiloMeter())

const AU = 1.49597870700e11 # m
const AU_ERR = 3 # m

"""
    astronomical_unit([[T=Float64], unit=KiloMeter()])

Return the astronomical unit ``au`` and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
astronomical_unit(::Type{Float64}, ::Meter) = AU, AU_ERR
astronomical_unit(::Type{Float64}, ::KiloMeter) = astronomical_unit(Float64, Meter()) .* 1e-3
astronomical_unit() = astronomical_unit(Float64, KiloMeter())

const LC = 1.48082686741e-8 # dimensionless
const LC_ERR = 2e-17 # dimensionless

"""
    drift_rate([T=Float64,] TCG, TCB)

Return the rate of drift between TCG and TCB ``L_C`` and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
drift_rate(
    ::Type{Float64},
    ::GeocentricCoordinateTime,
    ::BarycentricCoordinateTime,
) = LC, LC_ERR
drift_rate(
    ::GeocentricCoordinateTime,
    ::BarycentricCoordinateTime,
) = drift_rate(Float64, TCG, TCB)

const GM_S_TCB = 1.32712442099e20 # m³/s²
const GM_S_TCB_ERR = 1e10 # m³/s²
const GM_S_TDB = 1.32712440041e20 # m³/s²
const GM_S_TDB_ERR = 1e10 # m³/s²

"""
    heliocentric_gravitational_const([[[T=Float64], scale=TDB], unit=KiloMeter()])

Return the heliocentric gravitational constant ``GM_S`` and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
heliocentric_gravitational_const(
    ::Type{Float64},
    ::BarycentricCoordinateTime,
    ::Meter,
) = GM_S_TCB, GM_S_TCB_ERR
heliocentric_gravitational_const(
    ::Type{Float64},
    ::BarycentricDynamicalTime,
    ::Meter,
) = GM_S_TDB, GM_S_TDB_ERR
heliocentric_gravitational_const(
    ::Type{Float64},
    scale,
    ::KiloMeter,
) = heliocentric_gravitational_const(Float64, scale, Meter()) .* 1e-9
heliocentric_gravitational_const(s) = heliocentric_gravitational_const(Float64, s, KiloMeter())
heliocentric_gravitational_const() = heliocentric_gravitational_const(Float64, TDB, KiloMeter())

const AE = 6.3781366e6 # m
const AE_ERR = 1e-1 # m

"""
    equatorial_radius_earth([[T=Float64], unit=KiloMeter()])

Return the equatorial radius of Earth ``a\\_E`` and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
equatorial_radius_earth(::Type{Float64}, ::Meter) = AE, AE_ERR
equatorial_radius_earth(::Type{Float64}, ::KiloMeter) = equatorial_radius_earth(Float64, Meter()) .* 1e-3
equatorial_radius_earth() = equatorial_radius_earth(Float64, KiloMeter())

const J2 = 1.0826359e-3 # m
const J2_ERR = 1e-10 # m

"""
    j2_factor_earth([T=Float64])

Return ``J_2`` dynamical form factor of Earth and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
j2_factor_earth(::Type{Float64}) = J2, J2_ERR
j2_factor_earth() = j2_factor_earth(Float64)

const DJ2 = -3e-9 # 1/cy
const DJ2_ERR = 6e-10 # 1/cy

"""
    j2_rate_earth([T=Float64])

Return the time rate of change in ``J_2`` of Earth and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
j2_rate_earth(::Type{Float64}) = DJ2, DJ2_ERR
j2_rate_earth() = j2_rate_earth(Float64)

const GM_E_TCB = 3.986004418e14 # m³/s²
const GM_E_TCB_ERR = 8e5 # m³/s²
const GM_E_TT = 3.986004415e14 # m³/s²
const GM_E_TT_ERR = 8e5 # m³/s²
const GM_E_TDB = 3.986004356e14 # m³/s²
const GM_E_TDB_ERR = 8e5 # m³/s²

"""
    geocentric_gravitational_const([[[T=Float64], scale=TDB], unit=KiloMeter()])

Return the geocentric gravitational constant ``GM_E`` and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
geocentric_gravitational_const(
    ::Type{Float64},
    ::BarycentricCoordinateTime,
    ::Meter,
) = GM_E_TCB, GM_E_TCB_ERR
geocentric_gravitational_const(
    ::Type{Float64},
    ::TerrestrialTime,
    ::Meter,
) = GM_E_TT, GM_E_TT_ERR
geocentric_gravitational_const(
    ::Type{Float64},
    ::BarycentricDynamicalTime,
    ::Meter,
) = GM_E_TDB, GM_E_TDB_ERR
geocentric_gravitational_const(
    ::Type{Float64},
    scale,
    ::KiloMeter,
) = geocentric_gravitational_const(Float64, scale, Meter()) .* 1e-9
geocentric_gravitational_const(s) = geocentric_gravitational_const(Float64, s, KiloMeter())
geocentric_gravitational_const() = geocentric_gravitational_const(Float64, TDB, KiloMeter())

const W0 = 6.26368560e7 # m²/s²
const W0_ERR = 5e-1 # m²/s²

"""
    geoid_potential([T=Float64])

Return the geoid potential and the associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
geoid_potential(::Type{Float64}, ::Meter) = W0, W0_ERR
geoid_potential(::Type{Float64}, ::KiloMeter) = geoid_potential(Float64, Meter()) .* 1e-6
geoid_potential() = geoid_potential(Float64, KiloMeter())

const OMEGA = 7.292115e-5 # rad/s

"""
    mean_angular_velocity_earth([T=Float64])

Return nominal mean angular velocity of Earth ``\\omega``.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
mean_angular_velocity_earth(::Type{Float64}) = OMEGA
mean_angular_velocity_earth() = mean_angular_velocity_earth(Float64)

const MM_ME = 1.23000371e-2 # dimensionless
const MM_ME_ERR = 4e-10 # m²/s²

"""
    mass_ratio_moon([T=Float64])

Return the ratio of the mass of the Moon to the mass of Earth ``M_M/M_E`` and the
associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
mass_ratio_moon(::Type{Float64}) = MM_ME, MM_ME_ERR
mass_ratio_moon() = mass_ratio_moon(Float64)

const RATIOS = (
    ("mercury", "Me", 6.0236e6, 3e2),
    ("venus", "Ve", 4.08523719e5, 8e-3),
    ("mars", "Ma", 3.09870359e6, 2e-2),
    ("jupiter", "J", 1.047348644e3, 1.7e-5),
    ("saturn", "Sa", 3.4979018e3, 1e-4),
    ("uranus", "U", 2.290298e4, 3e-2),
    ("neptune", "N", 1.941226e4, 3e-2),
    ("pluto", "P", 1.36566e8, 2.8e4),
    ("eris", "Eris", 1.191e8, 1.4e6),
)

for (name, abbr, val, err) in RATIOS
    uname = uppercasefirst(name)
    fname = "mass_ratio_$name"
    func = Symbol(fname)
    @eval begin
        """
            $($fname)([T=Float64])

        Return the ratio of the mass of the Sun to the mass of $($uname) ``M_S/M_{$($abbr)}``
        and the associated uncertainty.

        # Reference

        - Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
            IAU working group on numerical standards for Fundamental Astronomy."
            *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
        """
        $func(::Type{Float64}) = $val, $err
        $func() = $func(Float64)

        export $func
    end
end

const RATIOS2 = (
    ("ceres", 4.72e-10, 3e-12),
    ("pallas", 1.03e-10, 3e-12),
    ("vesta", 1.35e-10, 3e-12),
)

for (name, val, err) in RATIOS2
    uname = uppercasefirst(name)
    fname = "mass_ratio_$name"
    func = Symbol(fname)
    @eval begin
        """
            $($fname)([T=Float64])

        Return the ratio of the mass of $($uname) to the mass of the Sun ``M_{$($uname)}/M_S``
        and the associated uncertainty.

        # Reference

        - Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
            IAU working group on numerical standards for Fundamental Astronomy."
            *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
        """
        $func(::Type{Float64}) = $val, $err
        $func() = $func(Float64)

        export $func
    end
end

const EPS = 8.4381406e4 # arcseconds
const EPS_ERR = 1e-3 # arcseconds

"""
    obliquity_j2000([T=Float64])

Return the ratio of the mass of the Moon to the mass of Earth ``M_M/M_E`` and the
associated uncertainty.

# Reference

- Luzum, Brian, et al. "The IAU 2009 system of astronomical constants: the report of the
    IAU working group on numerical standards for Fundamental Astronomy."
    *Celestial Mechanics and Dynamical Astronomy* 110.4 (2011): 293.
"""
obliquity_j2000(::Type{Float64}) = EPS, EPS_ERR
obliquity_j2000() = obliquity_j2000(Float64)

end

