using AstroTime: TDBEpoch, TTEpoch, centuries, j2000, value, julian_twopart
using ReferenceFrameRotations: angleaxis_to_dcm
using StaticArrays: SVector

using ..Bodies: CelestialBody,
    Mercury,
    Venus,
    EarthBarycenter,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    Sun,
    sun,
    grav_param
using ..Constants: astronomical_unit
using ..TwoBody: cartesian, transform, mean_anomaly, true_anomaly
using ..Util: sec2rad, normalize_angle

import ..Interfaces: AbstractEphemeris, position!, state!, velocity!

export SimonBretagnon, simon_bretagnon

struct SimonBretagnon <: AbstractEphemeris end
const simon_bretagnon = SimonBretagnon()

const ϵ = 0.40909280422232897

function semi_major(body::CelestialBody, t)
    throw(ArgumentError("Body '$body' is not supported by the Simon-Bretagnon ephemeris."))
end

semi_major(::Mercury, t) = 0.3870983098
eccentricity(::Mercury, t) = @evalpoly(t,
	0.2056317526,
	0.0002040653,
	-28349e-10,
)
inclination(::Mercury, t) = sec2rad(@evalpoly(t,
	7.00498625 * 3600,
	-214.25629,
	0.28977,
))
ascending_node(::Mercury, t) = sec2rad(@evalpoly(t,
	48.33089304 * 3600,
	-4515.21727,
	-31.79892,
))
longitude_of_perihelion(::Mercury, t) = sec2rad(@evalpoly(t,
	77.45611904 * 3600,
	5719.11590,
	-4.83016,
))
mean_longitude(::Mercury, t) = sec2rad(@evalpoly(t,
	252.25090552 * 3600,
	5381016286.88982,
	-1.92789,
))

semi_major(::Venus, t) = 0.7233298200
eccentricity(::Venus, t) = @evalpoly(t,
	0.0067719164,
	-0.0004776521,
	98127e-10,
)
inclination(::Venus, t) = sec2rad(@evalpoly(t,
	3.39466189 * 3600,
	-30.84437,
	-11.67836,
))
ascending_node(::Venus, t) = sec2rad(@evalpoly(t,
	76.67992019 * 3600,
	-10008.48154,
	-51.32614,
))
longitude_of_perihelion(::Venus, t) = sec2rad(@evalpoly(t,
	131.56370300 * 3600,
	175.48640,
	-498.48184,
))
mean_longitude(::Venus, t) = sec2rad(@evalpoly(t,
	181.97980085 * 3600,
	2106641364.33548,
	0.59381,
))

semi_major(::EarthBarycenter, t) = 1.0000010178
eccentricity(::EarthBarycenter, t) = @evalpoly(t,
	0.0167086342,
	-0.0004203654,
	-0.0000126734,
)
inclination(::EarthBarycenter, t) = sec2rad(@evalpoly(t,
	0.0,
	469.97289,
	-3.35053 ,
))
ascending_node(::EarthBarycenter, t) = sec2rad(@evalpoly(t,
	174.87317577 * 3600,
	-8679.27034,
	15.34191,
))
longitude_of_perihelion(::EarthBarycenter, t) = sec2rad(@evalpoly(t,
	102.93734808 * 3600,
	11612.35290,
	53.27577,
))
mean_longitude(::EarthBarycenter, t) = sec2rad(@evalpoly(t,
	100.46645683 * 3600,
	1295977422.83429,
	-2.04411,
))

semi_major(::Mars, t) = @evalpoly(t,
    1.5236793419,
    3e-10,
)
eccentricity(::Mars, t) = @evalpoly(t,
	0.0934006477,
	0.0009048438,
	-80641e-10,
)
inclination(::Mars, t) = sec2rad(@evalpoly(t,
	1.84972648 * 3600,
	-293.31722,
	-8.11830,
))
ascending_node(::Mars, t) = sec2rad(@evalpoly(t,
	49.55809321 * 3600,
	-10620.90088,
	-230.57416,
))
longitude_of_perihelion(::Mars, t) = sec2rad(@evalpoly(t,
	336.06023395 * 3600,
	15980.45908,
	-62.32800,
))
mean_longitude(::Mars, t) = sec2rad(@evalpoly(t,
	355.43299958 * 3600,
	689050774.93988,
	0.94264,
))

semi_major(::Jupiter, t) = @evalpoly(t,
    5.2026032092,
    19132e-10,
    -39e-10,
)
eccentricity(::Jupiter, t) = @evalpoly(t,
	0.0484979255,
	0.0016322542,
	-0.0000471366,
)
inclination(::Jupiter, t) = sec2rad(@evalpoly(t,
	1.30326698 * 3600,
	-71.55890,
	11.95297,
))
ascending_node(::Jupiter, t) = sec2rad(@evalpoly(t,
	100.46440702 * 3600,
	6362.03561,
	326.52178,
))
longitude_of_perihelion(::Jupiter, t) = sec2rad(@evalpoly(t,
	14.33120687 * 3600,
	7758.75163,
	259.95938,
))
mean_longitude(::Jupiter, t) = sec2rad(@evalpoly(t,
	34.35151874 * 3600,
	109256603.77991,
	-30.60378,
))

semi_major(::Saturn, t) = @evalpoly(t,
    9.5549091915,
    -0.0000213896,
    444e-10,
)
eccentricity(::Saturn, t) = @evalpoly(t,
	0.0555481426,
	-0.0034664062,
	-0.0000643639,
)
inclination(::Saturn, t) = sec2rad(@evalpoly(t,
	2.48887878 * 3600,
	91.85195,
	-17.66225,
))
ascending_node(::Saturn, t) = sec2rad(@evalpoly(t,
	113.66550252 * 3600,
	-9240.19942,
	-66.23743,
))
longitude_of_perihelion(::Saturn, t) = sec2rad(@evalpoly(t,
	93.05723748 * 3600,
	20395.49439,
	190.25952,
))
mean_longitude(::Saturn, t) = sec2rad(@evalpoly(t,
	50.07744430 * 3600,
	43996098.55732,
	75.61614,
))

semi_major(::Uranus, t) = @evalpoly(t,
    19.2184460618,
    -3716e-10,
    979e-10,
)
eccentricity(::Uranus, t) = @evalpoly(t,
	0.0463812221,
	-0.0002729293,
	0.0000078913,
)
inclination(::Uranus, t) = sec2rad(@evalpoly(t,
	0.77319689 * 3600,
	-60.72723,
	1.25759,
))
ascending_node(::Uranus, t) = sec2rad(@evalpoly(t,
	74.00595701 * 3600,
	2669.15033,
	145.93964,
))
longitude_of_perihelion(::Uranus, t) = sec2rad(@evalpoly(t,
	173.00529106 * 3600,
	3215.56238,
	-34.09288,
))
mean_longitude(::Uranus, t) = sec2rad(@evalpoly(t,
    314.05500511 * 3600,
    15424811.93933,
    -1.75083,
))

semi_major(::Neptune, t) = @evalpoly(t,
    30.1103868694,
    -16635e-10,
    686e-10,
)
eccentricity(::Neptune, t) = @evalpoly(t,
	0.0094557470,
	0.0000603263,
)
inclination(::Neptune, t) = sec2rad(@evalpoly(t,
	1.76995259 * 3600,
	8.12333,
	0.08135,
))
ascending_node(::Neptune, t) = sec2rad(@evalpoly(t,
	131.78405702 * 3600,
	-221.94322,
	-0.78728,
))
longitude_of_perihelion(::Neptune, t) = sec2rad(@evalpoly(t,
	48.12027554 * 3600,
	1050.71912,
	27.39717,
))
mean_longitude(::Neptune, t) = sec2rad(@evalpoly(t,
	304.34866548 * 3600,
	7865503.20744,
	0.21103,
))

# Tables for trigonometric terms to be added to the mean elements of
# the semi-major axes

kp(::Mercury) = SVector(69613, 75645, 88306, 59899, 15746, 71087, 142173, 3086, 0)
kp(::Venus) = SVector(21863, 32794, 26934, 10931, 26250, 43725, 53867, 28939, 0)
kp(::EarthBarycenter) = SVector(16002, 21863, 32004, 10931, 14529, 16368, 15318, 32794, 0)
kp(::Mars) = SVector(6345, 7818, 15636, 7077, 8184, 14163, 1107, 4872, 0)
kp(::Jupiter) = SVector(1760, 1454, 1167, 880, 287, 2640, 19, 2047, 1454)
kp(::Saturn) = SVector(574, 0, 880, 287, 19, 1760, 1167, 306, 574)
kp(::Uranus) = SVector(204, 0, 177, 1265, 4, 385, 200, 208, 204)
kp(::Neptune) = SVector(0, 102, 106, 4, 98, 1367, 487, 204, 0)

ca(::Mercury) = SVector(4, -13, 11, -9, -9, -3, -1, 4, 0)
ca(::Venus) = SVector(-156, 59, -42, 6, 19, -20, -10, -12, 0)
ca(::EarthBarycenter) = SVector(64, -152, 62, -8, 32, -41, 19, -11, 0)
ca(::Mars) = SVector(124, 621, -145, 208, 54, -57, 30, 15, 0)
ca(::Jupiter) = SVector(-23437, -2634, 6601, 6259, -1507, -1821, 2620, -2115, -1489)
ca(::Saturn) = SVector(62911, -119919, 79336, 17814, -24241, 12068, 8306, -4893, 8902)
ca(::Uranus) = SVector(389061, -262125, -44088, 8387, -22976, -2093, -615, -9720, 6633)
ca(::Neptune) = SVector(-412235, -157046, -31430, 37817, -9740, -13, -7449, 9644, 0)

sa(::Mercury) = SVector(-29, -1, 9, 6, -6, 5, 4, 0, 0)
sa(::Venus) = SVector(-48, -125, -26, -37, 18, -13, -20, -2, 0)
sa(::EarthBarycenter) = SVector(-150, -46, 68, 54, 14, 24, -28, 22, 0)
sa(::Mars) = SVector(-621, 532, -694, -20, 192, -94, 71, -73, 0)
sa(::Jupiter) = SVector(-14614, -19828, -5869, 1881, -4372, -2255, 782, 930, 913)
sa(::Saturn) = SVector(139737, 0, 24667, 51123, -5102, 7429, -4095, -1976, -9566)
sa(::Uranus) = SVector(-138081, 0, 37205, -49039, -41901, -33872, -27037, -12474, 18797)
sa(::Neptune) = SVector(0, 28492, 133236, 69654, 52322, -49577, -26430, -3593, 0)

# Tables giving the trigonometric terms to be added to the mean
# elements of the mean longitudes

kq(::Mercury) = SVector(3086, 15746, 69613, 59899, 75645, 88306, 12661, 2658, 0, 0)
kq(::Venus) = SVector(21863, 32794, 10931, 73, 4387, 26934, 1473, 2157, 0, 0)
kq(::EarthBarycenter) = SVector(10, 16002, 21863, 10931, 1473, 32004, 4387, 73, 0, 0)
kq(::Mars) = SVector(10, 6345, 7818, 1107, 15636, 7077, 8184, 532, 10, 0)
kq(::Jupiter) = SVector(19, 1760, 1454, 287, 1167, 880, 574, 2640, 19, 1454)
kq(::Saturn) = SVector(19, 574, 287, 306, 1760, 12, 31, 38, 19, 574)
kq(::Uranus) = SVector(4, 204, 177, 8, 31, 200, 1265, 102, 4, 204)
kq(::Neptune) = SVector(4, 102, 106, 8, 98, 1367, 487, 204, 4, 102)

cl(::Mercury) = SVector(21, -95, -157, 41, -5, 42, 23, 30, 0, 0)
cl(::Venus) = SVector(-160, -313, -235, 60, -74, -76, -27, 34, 0, 0)
cl(::EarthBarycenter) = SVector(-325, -322, -79, 232, -52, 97, 55, -41, 0, 0)
cl(::Mars) = SVector(2268, -979, 802, 602, -668, -33, 345, 201, -55, 0)
cl(::Jupiter) = SVector(7610, -4997, -7689, -5841, -2617, 1115, -748, -607, 6074, 354)
cl(::Saturn) = SVector(-18549, 30125, 20012, -730, 824, 23, 1289, -352, -14767, -2062)
cl(::Uranus) = SVector(-135245, -14594, 4197, -4030, -5630, -2898, 2540, -306, 2939, 1986)
cl(::Neptune) = SVector(89948, 2103, 8963, 2695, 3682, 1648, 866, -154, -1963, -283)

sl(::Mercury) = SVector(-342, 136, -23, 62, 66, -52, -33, 17, 0, 0)
sl(::Venus) = SVector(524, -149, -35, 117, 151, 122, -71, -62, 0, 0)
sl(::EarthBarycenter) = SVector(-105, -137, 258, 35, -116, -88, -112, -80, 0, 0)
sl(::Mars) = SVector(854, -205, -936, -240, 140, -341, -97, -232, 536, 0)
sl(::Jupiter) = SVector(-56980, 8016, 1012, 1448, -3024, -3710, 318, 503, 3767, 577)
sl(::Saturn) = SVector(138606, -13478, -4964, 1441, -1319, -1482, 427, 1236, -9167, -1918)
sl(::Uranus) = SVector(71234, -41116, 5334, -4935, -1848, 66, 434, -1748, 3780, -701)
sl(::Neptune) = SVector(-47645, 11647, 2166, 3194, 679, 0, -244, -419, -2531, 48)

function state!(pos, vel, ::SimonBretagnon, ep::Epoch, ::Sun, body::CelestialBody)
    t = value(centuries(j2000(TDBEpoch(ep)))) / 10.0

    a = semi_major(body, t)
    e = eccentricity(body, t)
    i = inclination(body, t)
    Ω = mod2pi(ascending_node(body, t))
    ω′ = mod2pi(longitude_of_perihelion(body, t))
    λm = mean_longitude(body, t)

    dμ = 0.35953620 * t;

    kp0 = kp(body)[1:8]
    kq0 = kq(body)[1:8]
    ca0 = ca(body)[1:8]
    sa0 = sa(body)[1:8]
    cl0 = cl(body)[1:8]
    sl0 = sl(body)[1:8]

    arga = kp0 .* dμ
    argl = kq0 .* dμ
    a += sum((ca0 .* cos.(arga) .+ sa0 .* sin.(arga)) * 1e-7)
    λm += sum((cl0 .* cos.(argl) .+ sl0 .* sin.(argl)) * 1e-7)

    kp1 = kp(body)[9:end]
    kq1 = kq(body)[9:end]
    ca1 = ca(body)[9:end]
    sa1 = sa(body)[9:end]
    cl1 = cl(body)[9:end]
    sl1 = sl(body)[9:end]

    arga = kp1 .* dμ
    argl = kq1 .* dμ
    a += sum(t * (ca1 .* cos.(arga) .+ sa1 .* sin.(arga)) * 1e-7)
    λm += sum(t * (cl1 .* cos.(argl) .+ sl1 .* sin.(argl)) * 1e-7)
    λm = normalize_angle(λm, 0.0)

    M = λm - ω′
    ω = ω′ - Ω
    ν = transform(mean_anomaly, true_anomaly, M, e)

    μ = grav_param(sun)

    r, v = cartesian(a * astronomical_unit(), e, i, Ω, ω, ν, μ)

    rot = angleaxis_to_dcm(-ϵ, [1, 0, 0])
    pos .+= rot * r
    vel .+= rot * v

    return pos, vel
end

function position!(pos, ::SimonBretagnon, ep::Epoch, ::Sun, body::CelestialBody)
    state!(pos, zeros(3), simon_bretagnon, ep, sun, body)
    return pos
end

function velocity!(vel, ::SimonBretagnon, ep::Epoch, ::Sun, body::CelestialBody)
    state!(zeros(3), vel, simon_bretagnon, ep, sun, body)
    return vel
end

function state!(pos, vel, ::SimonBretagnon, ep::Epoch, from::CelestialBody, to::CelestialBody)
    rv1 = state!(zeros(3), zeros(3), simon_bretagnon, ep, sun, from)
    rv2 = state!(zeros(3), zeros(3), simon_bretagnon, ep, sun, to)
    pos .+= rv2[1] .- rv1[1]
    vel .+= rv2[2] .- rv1[2]
    return pos, vel
end

function position!(pos, ::SimonBretagnon, ep::Epoch, from::CelestialBody, to::CelestialBody)
    state!(pos, zeros(3), simon_bretagnon, ep, from, to)
    return pos
end

function velocity!(vel, ::SimonBretagnon, ep::Epoch, from::CelestialBody, to::CelestialBody)
    state!(zeros(3), vel, simon_bretagnon, ep, from, to)
    return vel
end

