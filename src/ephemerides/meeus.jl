using AstroBase: AU, sec2rad
using AstroTime: TDBEpoch, TTEpoch, centuries, j2000, value, julian_twopart
using ReferenceFrameRotations: angleaxis_to_dcm
using ..Bodies: CelestialBody,
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    Sun,
    sun,
    grav_param
using ..EarthAttitude: obliquity_of_ecliptic_06
using ..TwoBody: cartesian, transform, mean_anomaly, true_anomaly

export Meeus, meeus

struct Meeus <: AbstractEphemeris end
const meeus = Meeus()

semi_major(::Mercury, t) = 0.387098310AU
eccentricity(::Mercury, t) = @evalpoly(t,
    0.20563175,
    0.000020406,
    -0.0000000284,
    -0.00000000017,
)
inclination(::Mercury, t) = deg2rad(@evalpoly(t,
    7.004986,
    -0.0059516,
    0.00000081,
    0.000000041,
))
ascending_node(::Mercury, t) = deg2rad(@evalpoly(t,
    48.330893,
    -0.1254229,
    -0.00008833,
    -0.000000196,
))
longitude_of_perihelion(::Mercury, t) = deg2rad(@evalpoly(t,
    77.456119,
    0.1588643,
    -0.00001343,
    0.000000039,
))
mean_longitude(::Mercury, t) = deg2rad(@evalpoly(t,
    252.250906,
    149472.6746358,
    -0.00000535,
    0.000000002,
))

semi_major(::Venus, t) = 0.723329820AU
eccentricity(::Venus, t) = @evalpoly(t,
    0.00677188,
    -0.000047766,
    0.0000000975,
    0.00000000044,
)
inclination(::Venus, t) = deg2rad(@evalpoly(t,
    3.394662,
    -0.0008568,
    -0.00003244,
    0.000000010,
))
ascending_node(::Venus, t) = deg2rad(@evalpoly(t,
    76.679920,
    -0.2780080,
    -0.00014256,
    -0.000000198,
))
longitude_of_perihelion(::Venus, t) = deg2rad(@evalpoly(t,
    131.563707,
    0.0048646,
    -0.00138232,
    -0.000005332,
))
mean_longitude(::Venus, t) = deg2rad(@evalpoly(t,
    181.979801,
    58517.8156760,
    0.00000165,
    -0.000000002,
))

semi_major(::Earth, t) = 1.000001018AU
eccentricity(::Earth, t) = @evalpoly(t,
    0.01670862,
    -0.000042037,
    -0.0000001236,
    0.00000000004,
)
inclination(::Earth, t) = deg2rad(@evalpoly(t,
    0.0,
    0.0130546,
    -0.00000931,
    -0.000000034,
))
ascending_node(::Earth, t) = deg2rad(@evalpoly(t,
    174.873174,
    -0.2410908,
    0.00004067,
    -0.000001327,
))
longitude_of_perihelion(::Earth, t) = deg2rad(@evalpoly(t,
    102.937348,
    0.3225557,
    0.00015026,
    0.000000478,
))
mean_longitude(::Earth, t) = deg2rad(@evalpoly(t,
    100.466449,
    35999.3728519,
    -0.00000568,
    0.0,
))

semi_major(::Mars, t) = 1.523679342AU
eccentricity(::Mars, t) = @evalpoly(t,
    0.09340062,
    0.000090483,
    -0.0000000806,
    -0.00000000035,
)
inclination(::Mars, t) = deg2rad(@evalpoly(t,
    1.849726,
    -0.0081479,
    -0.00002255,
    -0.000000027,
))
ascending_node(::Mars, t) = deg2rad(@evalpoly(t,
    49.558093,
    -0.2949846,
    -0.00063993,
    -0.000002143,
))
longitude_of_perihelion(::Mars, t) = deg2rad(@evalpoly(t,
    336.060234,
    0.4438898,
    -0.00017321,
    0.000000300,
))
mean_longitude(::Mars, t) = deg2rad(@evalpoly(t,
    355.433275,
    19140.2993313,
    0.00000261,
    -0.000000003,
))

semi_major(::Jupiter, t) = @evalpoly(t, 5.202603191, 0.0000001913) * AU
eccentricity(::Jupiter, t) = @evalpoly(t,
    0.04849485,
    0.000163244,
    -0.0000004719,
    -0.00000000197,
)
inclination(::Jupiter, t) = deg2rad(@evalpoly(t,
    1.303270,
    -0.0019872,
    0.00003318,
    0.000000092,
))
ascending_node(::Jupiter, t) = deg2rad(@evalpoly(t,
    100.464441,
    0.1766828,
    0.00090387,
    -0.000007032,
))
longitude_of_perihelion(::Jupiter, t) = deg2rad(@evalpoly(t,
    14.331309,
    0.2155525,
    0.00072252,
    -0.000004590,
))
mean_longitude(::Jupiter, t) = deg2rad(@evalpoly(t,
    34.351484,
    3034.9056746,
    -0.00008501,
    0.000000004,
))

semi_major(::Saturn, t) = @evalpoly(t, 9.554909596, -0.0000021389) * AU
eccentricity(::Saturn, t) = @evalpoly(t,
    0.05550862,
    -0.000346818,
    -0.0000006456,
    0.00000000338,
)
inclination(::Saturn, t) = deg2rad(@evalpoly(t,
    2.488878,
    0.0025515,
    -0.00004903,
    0.000000018,
))
ascending_node(::Saturn, t) = deg2rad(@evalpoly(t,
    113.665524,
    -0.2566649,
    -0.00018345,
    0.000000357,
))
longitude_of_perihelion(::Saturn, t) = deg2rad(@evalpoly(t,
    93.056787,
    0.5665496,
    0.00052809,
    0.000004882,
))
mean_longitude(::Saturn, t) = deg2rad(@evalpoly(t,
    50.077471,
    1222.1137943,
    0.00021004,
    -0.000000019,
))

semi_major(::Uranus, t) = @evalpoly(t, 19.218446062, -3.72e-11, 9.8e-12) * AU
eccentricity(::Uranus, t) = @evalpoly(t,
    0.04638122,
    -0.000027293,
    0.0000000789,
    2.447e-13,
)
inclination(::Uranus, t) = deg2rad(@evalpoly(t,
    0.773197,
    -0.0016869,
    0.00000349,
    0.000000016,
))
ascending_node(::Uranus, t) = deg2rad(@evalpoly(t,
    74.005957,
    0.0741431,
    0.00040539,
    0.000000119,
))
longitude_of_perihelion(::Uranus, t) = deg2rad(@evalpoly(t,
    173.005291,
    0.0893212,
    -0.00009470,
    0.000000414,
))
mean_longitude(::Uranus, t) = deg2rad(@evalpoly(t,
    314.055005,
    428.4669983,
    -0.00000486,
    5.99e-9,
))

semi_major(::Neptune, t) = @evalpoly(t, 30.110386869, -0.0000001663, 0.00000000069) * AU
eccentricity(::Neptune, t) = @evalpoly(t,
    0.00898809,
    0.000006408,
    -0.0000000008,
)
inclination(::Neptune, t) = deg2rad(@evalpoly(t,
    1.769952,
    0.0002257,
    0.00000023,
    -0.0,
))
ascending_node(::Neptune, t) = deg2rad(@evalpoly(t,
    131.784057,
    -0.0061651,
    -0.00000219,
    -0.000000078,
))
longitude_of_perihelion(::Neptune, t) = deg2rad(@evalpoly(t,
    48.123691,
    0.0291587,
    0.00007051,
    -0.0,
))
mean_longitude(::Neptune, t) = deg2rad(@evalpoly(t,
    304.348665,
    218.4862002,
    0.00000059,
    -0.000000002,
))

semi_major(::Pluto, t) = @evalpoly(t, 39.48168677, -0.00076912) * AU
eccentricity(::Pluto, t) = @evalpoly(t, 0.24880766, 0.00006465)
inclination(::Pluto, t) = @evalpoly(t, deg2rad(17.14175), sec2rad(11.07))
ascending_node(::Pluto, t) = @evalpoly(t, deg2rad(110.30347), sec2rad(-37.33))
longitude_of_perihelion(::Pluto, t) = @evalpoly(t, deg2rad(224.06676), sec2rad(-132.25))
mean_longitude(::Pluto, t) = @evalpoly(t, deg2rad(238.92881), sec2rad(522747.90))

function state!(pos, vel, ::Meeus, ep::TDBEpoch, ::Sun, body::CelestialBody)
    t = value(centuries(j2000(ep)))

    a = semi_major(body, t)
    e = eccentricity(body, t)
    i = inclination(body, t)
    Ω = ascending_node(body, t)
    ω′ = longitude_of_perihelion(body, t)
    λm = mean_longitude(body, t)

    M = λm - ω′
    ω = ω′ - Ω
    ν = transform(mean_anomaly, true_anomaly, M, e)

    μ = grav_param(sun)

    r, v = cartesian(a, e, i, Ω, ω, ν, μ)

    ϵ = obliquity_of_ecliptic_06(value.(julian_twopart(TTEpoch(ep)))...)
    rot = angleaxis_to_dcm(-ϵ, [1, 0, 0])
    pos .+= rot * r
    vel .+= rot * v

    pos, vel
end

function position!(pos, ::Meeus, ep::TDBEpoch, ::Sun, body::CelestialBody)
    state!(pos, zeros(3), meeus, ep, sun, body)
    pos
end

function velocity!(vel, ::Meeus, ep::TDBEpoch, ::Sun, body::CelestialBody)
    state!(zeros(3), vel, meeus, ep, sun, body)
    vel
end
