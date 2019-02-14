export MINOR_BODIES, MINOR_BODY_NAMES

const MINOR_BODY_NAMES = (
    :Pluto,
)

for minor in MINOR_BODY_NAMES
    @eval begin
        struct $minor <: MinorBody end
        export $minor
    end
end
@eval const MINOR_BODIES = [$(MINOR_BODY_NAMES...)]
