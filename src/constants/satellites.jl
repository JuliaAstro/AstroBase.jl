const SATELLITES = (
    :Moon,
)

for satellite in SATELLITES
    @eval begin
        struct $satellite <: NaturalSatellite end
        export $satellite
    end
end
