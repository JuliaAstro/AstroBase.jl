# AstroBase

*Interfaces, types, and functions for space science packages*

[![Build Status](https://github.com/JuliaAstro/AstroBase.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaAstro/AstroBase.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaAstro/AstroBase.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/AstroBase.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/AstroBase.jl/stable)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/AstroBase.jl/dev)

AstroBase.jl is a "thick" base package for developing space science libraries and solutions in Julia.
It provides fundamental algorithms and types as well as common interfaces that can be extended by
downstream packages.

# Features

AstroBase.jl provides a number submodules which can be individually imported,
e.g. `using AstroBase.Time`, and provide the following functionality:

- [Astronomical Time (`Time`)](https://juliaastro.github.io/AstroBase.jl/dev/modules/time.html):
	A wrapper for the [AstroTime.jl](https://github.com/JuliaAstro/AstroTime.jl) package

