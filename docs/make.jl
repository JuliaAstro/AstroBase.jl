using Documenter, AstroBase

setup = quote
    using AstroBase
end
DocMeta.setdocmeta!(AstroBase, :DocTestSetup, setup; recursive = true)

makedocs(
    format = Documenter.HTML(
	prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [AstroBase],
    sitename = "AstroBase.jl",
    authors = "Helge Eichhorn and the AstroBase.jl contributors",
    pages = [
	"Home" => "index.md",
	"Modules" => [
	    "Celestial Bodies" => joinpath("modules", "bodies.md"),
	    "Constants" => joinpath("modules", "constants.md"),
	    "Coordinates" => joinpath("modules", "coords.md"),
	    "Earth Attitude" => joinpath("modules", "earth_attitude.md"),
	    "Ephemerides" => joinpath("modules", "ephemerides.md"),
	    "Reference Frames" => joinpath("modules", "frames.md"),
	    "Time" => joinpath("modules", "time.md"),
	    "Two-Body Problem" => joinpath("modules", "two_body.md"),
	    "Utilities" => joinpath("modules", "util.md"),
	],
    ],
    # strict = true,
    doctest = false,
)

deploydocs(
    repo = "github.com/JuliaAstro/AstroBase.jl.git",
    push_preview = true,
)
