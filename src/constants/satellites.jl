const SATELLITES = (
    :Moon,
)

for satellite in SATELLITES
    @eval begin
        struct $satellite <: NaturalSatellite end
        export $satellite
    end
end

μ(::Type{Moon}) = 4902.80007km^3/s^2
j2(::Type{Moon}) = 203.21568e-6
mean_radius(::Type{Moon}) = 1737.4km
subplanetary_equatorial_radius(::Type{Moon}) = 1737.4km
along_orbit_equatorial_radius(::Type{Moon}) = 1737.4km
polar_radius(::Type{Moon}) = 1737.4km
deviation(::Type{Moon}) = 2.5km
maximum_elevation(::Type{Moon}) = 7.5km
maximum_depression(::Type{Moon}) = 5.6km
naif_id(::Type{Moon}) = 301
parent(::Type{Moon}) = EarthBarycenter
ra₀(::Type{Moon}) = deg2rad(269.9949)
ra₁(::Type{Moon}) = deg2rad(0.0031)
dec₀(::Type{Moon}) = deg2rad(66.5392)
dec₁(::Type{Moon}) = deg2rad(0.0130)
w₀(::Type{Moon}) = deg2rad(38.3213)
w₁(::Type{Moon}) = deg2rad(13.17635815)
w₂(::Type{Moon}) = deg2rad(-1.4e-12)
a(::Type{Moon}) = deg2rad.([-3.8787, -0.1204, 0.07, -0.0172, 0.0, 0.0072, 0.0, 0.0, 0.0, -0.0052, 0.0, 0.0, 0.0043])
d(::Type{Moon}) = deg2rad.([1.5419, 0.0239, -0.0278, 0.0068, 0.0, -0.0029, 0.0009, 0.0, 0.0, 0.0008, 0.0, 0.0, -0.0009])
w(::Type{Moon}) = deg2rad.([3.5610, 0.1208, -0.0642, 0.0158, 0.0252, -0.0066, -0.0047, -0.0046, 0.0028, 0.0052, 0.004, 0.0019, -0.0044])
θ₀(::Type{Moon}) = deg2rad.([125.045, 250.089, 260.008, 176.625, 357.529, 311.589, 134.963, 276.617, 34.226, 15.134, 119.743, 239.961, 25.053])
θ₁(::Type{Moon}) = deg2rad.([-0.0529921, -0.1059842, 13.0120009, 13.3407154, 0.9856003, 26.4057084, 13.0649930, 0.3287146, 1.7484877, -0.1589763, 0.0036096, 0.1643573, 12.9590088])*JULIAN_CENTURY
