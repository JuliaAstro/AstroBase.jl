# AstroBase

*Interfaces, types, and functions for space science packages*

[![Build Status](https://github.com/JuliaAstro/AstroBase.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaAstro/AstroBase.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaAstro/AstroBase.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/AstroBase.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/AstroBase.jl/stable)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/AstroBase.jl/dev)

AstroBase.jl is a "thick" base package for developing space science libraries and solutions in Julia.
It provides fundamental algorithms and types as well as common interfaces that can be extended by
downstream packages.

## Features

AstroBase.jl provides a number of submodules which can be individually imported,
e.g. `using AstroBase.Time`, and provide the following functionality:

<!-- FIXME: Use link to stable docs once this is about to be published -->
- [`Bodies`](https://juliaastro.github.io/AstroBase.jl/dev/modules/bodies/):
	Types representing celestial bodies and associated constants
- [`Constants`](https://juliaastro.github.io/AstroBase.jl/dev/modules/constants/):
	Astronomical constants
- [`Coords`](https://juliaastro.github.io/AstroBase.jl/dev/modules/coords/):
	Coordinate representations of objects in space (Cartesian, Keplerian) and transformations
- [`EarthAttitude`](https://juliaastro.github.io/AstroBase.jl/dev/modules/earth_attitude/):
	Earth attitude modelling tools
- [`Ephemerides`](https://juliaastro.github.io/AstroBase.jl/dev/modules/ephemerides/):
	Semi-analytical planetary ephemerides
- [`Frames`](https://juliaastro.github.io/AstroBase.jl/dev/modules/frames/):
	Types representing quasi-inertial and rotating references frames and the associates
	transformations
- [`Interfaces`](https://juliaastro.github.io/AstroBase.jl/dev/modules/interfaces/):
	Abstract types and base functions
- [`Time`](https://juliaastro.github.io/AstroBase.jl/dev/modules/time/):
	A wrapper for the [AstroTime.jl](https://github.com/JuliaAstro/AstroTime.jl) package
- [`TwoBody`](https://juliaastro.github.io/AstroBase.jl/dev/modules/two_body/):
	Functions related to the two-body problem of celestial mechanics
- [`Util`](https://juliaastro.github.io/AstroBase.jl/dev/modules/util/):
	Various utility functions e.g. angle conversions

## Packages Using AstroBase.jl

- [JPLEphemeris.jl](https://github.com/JuliaAstro/JPLEphemeris.jl)

## Documentation

Please refer to the [documentation](https://JuliaAstro.github.io/AstroBase.jl/stable/)
for additional information

