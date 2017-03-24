export MINOR_BODIES

const MINOR_BODIES = [
    :Pluto,
]

for minor in MINOR_BODIES
    @eval begin
        struct $minor <: MinorBody end
        export $minor
    end
end
