# Celestial Bodies

```@meta
DocTestSetup = quote
    using AstroBase
end
```

```@contents
Pages = ["bodies.md"]
```

## Functions

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Order = [:function]
```

## Macros

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Order = [:macro]
```

## Base Types

```@docs
NAIFId
CelestialBody
Barycenter
Planet
MinorBody
NaturalSatellite
EarthSatellite
MarsSatellite
JupiterSatellite
SaturnSatellite
UranusSatellite
NeptuneSatellite
PlutoSatellite
```

## Sun and Solar System Barycenter

```@docs
sun
Sun
ssb
SolarSystemBarycenter
```

## Planets and Planetary Barycenters

```@docs
mercury
Mercury
mercury_barycenter
MercuryBarycenter
venus
Venus
venus_barycenter
VenusBarycenter
earth
Earth
earth_barycenter
EarthBarycenter
mars
Mars
mars_barycenter
MarsBarycenter
jupiter
Jupiter
jupiter_barycenter
JupiterBarycenter
saturn
Saturn
saturn_barycenter
SaturnBarycenter
uranus
Uranus
uranus_barycenter
UranusBarycenter
neptune
Neptune
neptune_barycenter
NeptuneBarycenter
```

## Minor Bodies

```@docs
pluto_barycenter
PlutoBarycenter
```

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: MinorBody) || t isa MinorBody
```

## Natural Satellites

### Earth's Moon

```@docs
luna
moon
Luna
```

### Mars Satellites

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: MarsSatellite) || t isa MarsSatellite
```

### Jupiter Satellites

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: JupiterSatellite) || t isa JupiterSatellite
```

### Saturn Satellites

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: SaturnSatellite) || t isa SaturnSatellite
```

### Uranus Satellites

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: UranusSatellite) || t isa UranusSatellite
```

### Neptune Satellites

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: NeptuneSatellite) || t isa NeptuneSatellite
```

### Pluto Satellites

```@autodocs
Modules = [AstroBase.Bodies]
Private = false
Filter = t -> (isconcretetype(t) && t <: PlutoSatellite) || t isa PlutoSatellite
```

```@meta
DocTestSetup = nothing
```

