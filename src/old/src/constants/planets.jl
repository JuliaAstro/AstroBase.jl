# Luzum, Brian, et al. "The IAU 2009 system of astronomical constants:
# the report of the IAU working group on numerical standards for Fundamental
# Astronomy." Celestial Mechanics and Dynamical Astronomy 110.4 (2011): 293-304.
# Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic
# coordinates and rotational elements: 2009." Celestial Mechanics and Dynamical
# Astronomy 109.2 (2011): 101-135.
export PLANETS, PLANET_NAMES

const PLANET_NAMES = (
    :Mercury,
    :Venus,
    :Earth,
    :Mars,
    :Jupiter,
    :Saturn,
    :Uranus,
    :Neptune
)

for (i, planet) in enumerate(PLANET_NAMES)
    barycenter = Symbol(planet, "Barycenter")
    @eval begin
        struct $planet <: Planet end
        struct $barycenter <: Barycenter end
        export $planet, $barycenter

        naif_id(::Type{$barycenter}) = $i
        parent(::Type{$barycenter}) = SSB
        body(::Type{$barycenter}) = $planet
    end
end
@eval const PLANETS = [$(PLANET_NAMES...)]

μ(::Type{Mercury}) = 2.2032e4
j2(::Type{Mercury}) = 2.027e-4
mean_radius(::Type{Mercury}) = 2439.7
equatorial_radius(::Type{Mercury}) = 2439.7
polar_radius(::Type{Mercury}) = 2439.7
deviation(::Type{Mercury}) = 1.0
maximum_elevation(::Type{Mercury}) = 4.6
maximum_depression(::Type{Mercury}) = 2.5
naif_id(::Type{Mercury}) = 199
parent(::Type{Mercury}) = MercuryBarycenter
α₀(::Type{Mercury}) = deg2rad(281.0097)
α₁(::Type{Mercury}) = deg2rad(-0.0328)
δ₀(::Type{Mercury}) = deg2rad(61.41430)
δ₁(::Type{Mercury}) = deg2rad(-0.0049)
w₀(::Type{Mercury}) = deg2rad(329.5469)
w₁(::Type{Mercury}) = deg2rad(6.1385025)
w(::Type{Mercury}) = deg2rad.([0.00993822, -0.00104581, -0.00010280, -0.00002364,
    -0.00000532])
θ₀(::Type{Mercury}) = deg2rad.([174.791086, 349.582171, 164.373257, 339.164343,
    153.955429])
θ₁(::Type{Mercury}) = deg2rad.([4.092335, 8.184670, 12.277005, 16.369340,
    20.461675])*36525.0

μ(::Type{Venus}) = 3.24859e5
j2(::Type{Venus}) = 6e-5
mean_radius(::Type{Venus}) = 6051.8
equatorial_radius(::Type{Venus}) = 6051.8
polar_radius(::Type{Venus}) = 6051.8
deviation(::Type{Venus}) = 1.0
maximum_elevation(::Type{Venus}) = 11.0
maximum_depression(::Type{Venus}) = 2.0
naif_id(::Type{Venus}) = 299
parent(::Type{Venus}) = VenusBarycenter
α₀(::Type{Venus}) = deg2rad(272.76)
δ₀(::Type{Venus}) = deg2rad(67.16)
w₀(::Type{Venus}) = deg2rad(160.20)
w₁(::Type{Venus}) = deg2rad(-1.4813688)

μ(::Type{Earth}) = 3.986004418e5
j2(::Type{Earth}) = 1.08262668e-3
mean_radius(::Type{Earth}) = 6371.0084
equatorial_radius(::Type{Earth}) = 6378.1366
polar_radius(::Type{Earth}) = 6356.7519
deviation(::Type{Earth}) = 3.57
maximum_elevation(::Type{Earth}) = 8.85
maximum_depression(::Type{Earth}) = 11.52
naif_id(::Type{Earth}) = 399
parent(::Type{Earth}) = EarthBarycenter
α₁(::Type{Earth}) = deg2rad(-0.641)
δ₀(::Type{Earth}) = deg2rad(90)
δ₁(::Type{Earth}) = deg2rad(-0.557)
w₀(::Type{Earth}) = deg2rad(190.147)
w₁(::Type{Earth}) = deg2rad(360.9856235)

μ(::Type{Mars}) = 4.282837e4
j2(::Type{Mars}) = 1.964e-3
mean_radius(::Type{Mars}) = 3389.5
equatorial_radius(::Type{Mars}) = 3396.19
polar_radius(::Type{Mars}) = 3376.20
deviation(::Type{Mars}) = 3.0
maximum_elevation(::Type{Mars}) = 22.64
maximum_depression(::Type{Mars}) = 7.55
naif_id(::Type{Mars}) = 499
parent(::Type{Mars}) = MarsBarycenter
α₀(::Type{Mars}) = deg2rad(317.68143)
α₁(::Type{Mars}) = deg2rad(-0.1061)
δ₀(::Type{Mars}) = deg2rad(52.886500)
δ₁(::Type{Mars}) = deg2rad(-0.0609)
w₀(::Type{Mars}) = deg2rad(176.630)
w₁(::Type{Mars}) = deg2rad(350.89198226)

μ(::Type{Jupiter}) = 1.26686534e8
j2(::Type{Jupiter}) = 1.475e-2
mean_radius(::Type{Jupiter}) = 69911.0
equatorial_radius(::Type{Jupiter}) = 71492.0
polar_radius(::Type{Jupiter}) = 66854.0
deviation(::Type{Jupiter}) = 62.1
maximum_elevation(::Type{Jupiter}) = 31.0
maximum_depression(::Type{Jupiter}) = 102.0
naif_id(::Type{Jupiter}) = 599
parent(::Type{Jupiter}) = JupiterBarycenter
α₀(::Type{Jupiter}) = deg2rad(268.056595)
α₁(::Type{Jupiter}) = deg2rad(-0.006499)
δ₀(::Type{Jupiter}) = deg2rad(64.4953030)
δ₁(::Type{Jupiter}) = deg2rad(0.002413)
w₀(::Type{Jupiter}) = deg2rad(284.95)
w₁(::Type{Jupiter}) = deg2rad(870.536)
α(::Type{Jupiter}) = deg2rad.([0.000117, 0.000938, 0.001432, 0.00003, 0.002150])
δ(::Type{Jupiter}) = deg2rad.([0.00005, 0.000404, 0.000617, -0.000013, 0.000926])
θ₀(::Type{Jupiter}) = deg2rad.([99.360714, 175.895369, 300.323162, 114.012305,
    49.511251])
θ₁(::Type{Jupiter}) = deg2rad.([4850.4046, 1191.9605, 262.5475, 6070.2476, 64.3])

μ(::Type{Saturn}) = 3.7931187e7
j2(::Type{Saturn}) = 1.645e-2
mean_radius(::Type{Saturn}) = 58232.0
equatorial_radius(::Type{Saturn}) = 60268.0
polar_radius(::Type{Saturn}) = 54364.0
deviation(::Type{Saturn}) = 102.9
maximum_elevation(::Type{Saturn}) = 8.0
maximum_depression(::Type{Saturn}) = 205.0
naif_id(::Type{Saturn}) = 699
parent(::Type{Saturn}) = SaturnBarycenter
α₀(::Type{Saturn}) = deg2rad(40.589)
α₁(::Type{Saturn}) = deg2rad(-0.036)
δ₀(::Type{Saturn}) = deg2rad(83.537)
δ₁(::Type{Saturn}) = deg2rad(-0.004)
w₀(::Type{Saturn}) = deg2rad(38.90)
w₁(::Type{Saturn}) = deg2rad(810.7939024)

μ(::Type{Uranus}) = 5.793939e6
j2(::Type{Uranus}) = 1.2e-2
mean_radius(::Type{Uranus}) = 25362.0
equatorial_radius(::Type{Uranus}) = 25559.0
polar_radius(::Type{Uranus}) = 24973.0
deviation(::Type{Uranus}) = 16.8
maximum_elevation(::Type{Uranus}) = 28.0
maximum_depression(::Type{Uranus}) = 0.0
naif_id(::Type{Uranus}) = 799
parent(::Type{Uranus}) = UranusBarycenter
α₀(::Type{Uranus}) = deg2rad(257.311)
δ₀(::Type{Uranus}) = deg2rad(-15.1750)
w₀(::Type{Uranus}) = deg2rad(203.81)
w₁(::Type{Uranus}) = deg2rad(-501.1600928)

μ(::Type{Neptune}) = 6.836529e6
j2(::Type{Neptune}) = 4e-3
mean_radius(::Type{Neptune}) = 24622.0
equatorial_radius(::Type{Neptune}) = 24764.0
polar_radius(::Type{Neptune}) = 24341.0
deviation(::Type{Neptune}) = 8.0
maximum_elevation(::Type{Neptune}) = 14.0
maximum_depression(::Type{Neptune}) = 0.0
naif_id(::Type{Neptune}) = 899
parent(::Type{Neptune}) = NeptuneBarycenter
α₀(::Type{Neptune}) = deg2rad(299.36)
δ₀(::Type{Neptune}) = deg2rad(43.460)
w₀(::Type{Neptune}) = deg2rad(253.18)
w₁(::Type{Neptune}) = deg2rad(536.3128492)
α(::Type{Neptune}) = deg2rad.([0.7])
δ(::Type{Neptune}) = deg2rad.([-0.51])
w(::Type{Neptune}) = deg2rad.([-0.48])
θ₀(::Type{Neptune}) = deg2rad.([357.85])
θ₁(::Type{Neptune}) = deg2rad.([52.316])
