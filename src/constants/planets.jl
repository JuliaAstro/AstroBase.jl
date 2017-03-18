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

for planet in PLANETS
    @eval begin
        struct $planet <: Planet end
        export $planet
    end
end

μ(::Type{Earth}) = 3.986004418e5km^3/s^2
