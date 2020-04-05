var documenterSearchIndex = {"docs":
[{"location":"modules/constants.html#Constants-1","page":"Constants","title":"Constants","text":"","category":"section"},{"location":"modules/constants.html#","page":"Constants","title":"Constants","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/constants.html#","page":"Constants","title":"Constants","text":"Modules = [AstroBase.Constants]\nPrivate = false","category":"page"},{"location":"modules/constants.html#","page":"Constants","title":"Constants","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/coords.html#Coordinates-1","page":"Coordinates","title":"Coordinates","text":"","category":"section"},{"location":"modules/coords.html#","page":"Coordinates","title":"Coordinates","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/coords.html#","page":"Coordinates","title":"Coordinates","text":"Modules = [AstroBase.Coords]\nPrivate = false","category":"page"},{"location":"modules/coords.html#","page":"Coordinates","title":"Coordinates","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/interfaces.html#Interfaces-1","page":"Interfaces","title":"Interfaces","text":"","category":"section"},{"location":"modules/interfaces.html#","page":"Interfaces","title":"Interfaces","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/interfaces.html#","page":"Interfaces","title":"Interfaces","text":"Modules = [AstroBase.Interfaces]\nPrivate = false","category":"page"},{"location":"modules/interfaces.html#","page":"Interfaces","title":"Interfaces","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/earth_attitude.html#Earth-Attitude-1","page":"Earth Attitude","title":"Earth Attitude","text":"","category":"section"},{"location":"modules/earth_attitude.html#","page":"Earth Attitude","title":"Earth Attitude","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/earth_attitude.html#","page":"Earth Attitude","title":"Earth Attitude","text":"Modules = [AstroBase.EarthAttitude]\nPrivate = false","category":"page"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau1980","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau1980","text":"`iau1980`\n\nThe singleton instance of type IAU1980, representing the IAU 1980 family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau1982","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau1982","text":"`iau1982`\n\nThe singleton instance of type IAU1982, representing the IAU 1982 family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau1994","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau1994","text":"`iau1994`\n\nThe singleton instance of type IAU1994, representing the IAU 1994 family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau2000","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau2000","text":"`iau2000`\n\nThe singleton instance of type IAU2000, representing the IAU 2000 family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau2000a","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau2000a","text":"`iau2000a`\n\nThe singleton instance of type IAU2000A, representing the IAU 2000A family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau2000b","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau2000b","text":"`iau2000b`\n\nThe singleton instance of type IAU2000B, representing the IAU 2000B family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau2006","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau2006","text":"`iau2006`\n\nThe singleton instance of type IAU2006, representing the IAU 2006 family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.iau2006a","page":"Earth Attitude","title":"AstroBase.EarthAttitude.iau2006a","text":"`iau2006a`\n\nThe singleton instance of type IAU2006a, representing the IAU 2006A family of models.\n\n\n\n\n\n","category":"constant"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.cip_coords-Tuple{AbstractArray{T,2} where T}","page":"Earth Attitude","title":"AstroBase.EarthAttitude.cip_coords","text":"cip_coords(rbpn)\n\nExtract from the bias-precession-nutation matrix the X,Y coordinates of the Celestial Intermediate Pole.\n\nReferences\n\nERFA\n\n\n\n\n\n","category":"method"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.nutation","page":"Earth Attitude","title":"AstroBase.EarthAttitude.nutation","text":"nutation(model, ep)\n\nReturn the nutation components for a given epoch and model.\n\nArguments\n\nmodel: IAU model, one of: iau1980, iau2000a, iau2000b, iau2006\nep: An epoch\n\nOutput\n\nδψ: Nutation in longitude (radians)\nδϵ: Nutation in obliquity (radians)\n\nExample\n\njulia> ep = TTEpoch(2020, 1, 1)\n2020-01-01T00:00:00.000 TT\n\njulia> nutation(iau2006, ep)\n(-7.996558232098883e-5, -8.25141288270117e-6)\n\nReferences\n\nSOFA\n\n\n\n\n\n","category":"function"},{"location":"modules/earth_attitude.html#AstroBase.EarthAttitude.obliquity","page":"Earth Attitude","title":"AstroBase.EarthAttitude.obliquity","text":"obliquity(model, ep)\n\nReturn the mean obliquity of the ecliptic for a given epoch and model.\n\nArguments\n\nmodel: IAU model, one of: iau1980, iau2006\nep: An epoch\n\nOutput\n\nReturns the angle between the ecliptic and mean equator of date in radians.\n\nExample\n\njulia> ep = TTEpoch(2020, 1, 1)\n2020-01-01T00:00:00.000 TT\n\njulia> obliquity(iau2006, ep)\n0.40904718953841473\n\nReferences\n\nSOFA\n\n\n\n\n\n","category":"function"},{"location":"modules/earth_attitude.html#","page":"Earth Attitude","title":"Earth Attitude","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/ephemerides.html#Ephemerides-1","page":"Ephemerides","title":"Ephemerides","text":"","category":"section"},{"location":"modules/ephemerides.html#","page":"Ephemerides","title":"Ephemerides","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/ephemerides.html#","page":"Ephemerides","title":"Ephemerides","text":"Modules = [AstroBase.Ephemerides]\nPrivate = false","category":"page"},{"location":"modules/ephemerides.html#","page":"Ephemerides","title":"Ephemerides","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/two_body.html#Two-Body-Problem-1","page":"Two-Body Problem","title":"Two-Body Problem","text":"","category":"section"},{"location":"modules/two_body.html#","page":"Two-Body Problem","title":"Two-Body Problem","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/two_body.html#","page":"Two-Body Problem","title":"Two-Body Problem","text":"Modules = [AstroBase.TwoBody]\nPrivate = false","category":"page"},{"location":"modules/two_body.html#AstroBase.TwoBody.transform-NTuple{4,Any}","page":"Two-Body Problem","title":"AstroBase.TwoBody.transform","text":"transform(from, to, a, ecc)\n\nTransform anomaly a from one anomaly type to another for an orbit with eccentricity ecc.\n\nArguments\n\nfrom, to: Anomaly types\ntrue_anomaly\neccentric_anomaly\nmean_anomaly\na: Current value of the anomaly\necc: Eccentricity of the orbit\n\nOutput\n\nReturns the transformed anomaly.\n\nReferences\n\nFarnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.   \"Robust resolution of Kepler’s equation in all eccentricity regimes.\"   Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.\n\n\n\n\n\n","category":"method"},{"location":"modules/two_body.html#","page":"Two-Body Problem","title":"Two-Body Problem","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/util.html#Utilities-1","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"modules/util.html#","page":"Utilities","title":"Utilities","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/util.html#","page":"Utilities","title":"Utilities","text":"Modules = [AstroBase.Util]\nPrivate = false","category":"page"},{"location":"modules/util.html#AstroBase.Util.rad2sec-Tuple{Any}","page":"Utilities","title":"AstroBase.Util.rad2sec","text":"rad2sec(rad)\n\nConvert an angle in radians to arcseconds.\n\nExample\n\njulia> rad2sec(0.5235987755982988)\n107999.99999999999\n\n\n\n\n\n","category":"method"},{"location":"modules/util.html#AstroBase.Util.sec2rad-Tuple{Any}","page":"Utilities","title":"AstroBase.Util.sec2rad","text":"sec2rad(sec)\n\nConvert an angle in arcseconds to radians.\n\nExample\n\njulia> sec2rad(3600 * 30)\n0.5235987755982988\n\n\n\n\n\n","category":"method"},{"location":"modules/util.html#AstroBase.Util.spherical_to_cartesian-Tuple{Any,Any}","page":"Utilities","title":"AstroBase.Util.spherical_to_cartesian","text":"sphericaltocartesian(theta, phi)\n\nConvert spherical coordinates to Cartesian.\n\nArguments\n\ntheta: longitude angle in radians\nphi: latitude angle in radians\n\nReturns\n\nx: magnitude of projection on x axis\ny: magnitude of projection on y axis\nz: magnitude of projection on z axis\n\nReferences\n\nERFA\n\n\n\n\n\n","category":"method"},{"location":"modules/util.html#","page":"Utilities","title":"Utilities","text":"DocTestSetup = nothing","category":"page"},{"location":"modules/frames.html#Reference-Frames-1","page":"Reference Frames","title":"Reference Frames","text":"","category":"section"},{"location":"modules/frames.html#","page":"Reference Frames","title":"Reference Frames","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/frames.html#","page":"Reference Frames","title":"Reference Frames","text":"Modules = [AstroBase.Frames]\nPrivate = false","category":"page"},{"location":"modules/frames.html#","page":"Reference Frames","title":"Reference Frames","text":"DocTestSetup = nothing","category":"page"},{"location":"index.html#AstroBase.jl-1","page":"Home","title":"AstroBase.jl","text":"","category":"section"},{"location":"modules/bodies.html#Celestial-Bodies-1","page":"Celestial Bodies","title":"Celestial Bodies","text":"","category":"section"},{"location":"modules/bodies.html#","page":"Celestial Bodies","title":"Celestial Bodies","text":"DocTestSetup = quote\n    using AstroBase\nend","category":"page"},{"location":"modules/bodies.html#","page":"Celestial Bodies","title":"Celestial Bodies","text":"Modules = [AstroBase.Bodies]\nPrivate = false","category":"page"},{"location":"modules/bodies.html#","page":"Celestial Bodies","title":"Celestial Bodies","text":"DocTestSetup = nothing","category":"page"}]
}