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
    ],
    # strict = true,
    doctest = false,
)

deploydocs(
    repo = "github.com/JuliaAstro/AstroBase.jl.git",
    push_preview = true,
)
