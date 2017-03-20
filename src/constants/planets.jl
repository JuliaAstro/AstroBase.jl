const PLANETS = (
    :Mercury,
    :Venus,
    :Earth,
    :Mars,
    :Jupiter,
    :Saturn,
    :Uranus,
    :Neptune
)

for (i, planet) in enumerate(PLANETS)
    barycenter = Symbol(planet, "Barycenter")
    @eval begin
        struct $planet <: Planet end
        struct $barycenter <: Barycenter end
        export $planet, $barycenter

        naif_id(::Type{$barycenter}) = $i
        parent(::Type{$barycenter}) = SSB
    end
end

μ(::Type{Earth}) = 3.986004418e5km^3/s^2
j2(::Type{Earth}) = 1.08262668e-3
mean_radius(::Type{Earth}) = 6371.0084km
equatorial_radius(::Type{Earth}) = 6378.1366km
polar_radius(::Type{Earth}) = 6356.7519km
deviation(::Type{Earth}) = 3.57km
maximum_elevation(::Type{Earth}) = 8.85km
maximum_depression(::Type{Earth}) = 11.52km
naif_id(::Type{Earth}) = 399
parent(::Type{Earth}) = EarthBarycenter
ra₀(::Type{Earth}) = 0.0
ra₁(::Type{Earth}) = deg2rad(-0.641)
ra₂(::Type{Earth}) = 0.0
dec₀(::Type{Earth}) = deg2rad(90)
dec₁(::Type{Earth}) = deg2rad(-0.557)
dec₂(::Type{Earth}) = 0.0
w₀(::Type{Earth}) = deg2rad(190.147)
w₁(::Type{Earth}) = deg2rad(360.9856235)
w₂(::Type{Earth}) = 0.0
a(::Type{Earth}) = [0.0]
d(::Type{Earth}) = [0.0]
w(::Type{Earth}) = [0.0]
θ₀(::Type{Earth}) = [0.0]
θ₁(::Type{Earth}) = [0.0]
