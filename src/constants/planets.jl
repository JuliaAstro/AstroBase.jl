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
