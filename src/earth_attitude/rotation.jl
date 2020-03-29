using AstroTime: SECONDS_PER_DAY, TDB, TT, UT1Epoch, fractionofday, j2000, value, julian_period, UT1
using ..Util: normalize_angle

export apparent_sidereal, earth_rotation_angle, equinoxes, mean_sidereal

function earth_rotation_angle(::IAU2000, ep)
    t = julian_period(ep; scale=UT1, raw=true)
    f = t % 1.0
    return mod2pi(2π * (f + 0.7790572732640 + 0.00273781191135448t))
end

include(joinpath("constants", "equinoxes.jl"))

function equinoxes(::IAU1994, ep::Epoch; scale=TDB)
    t = julian_period(ep; scale=scale, unit=centuries, raw=true)

    ω₀ = sec2rad(@evalpoly(t, 450160.280, -482890.539, 7.455, 0.008))
    ω = mod2pi(ω₀ + (-5.0t % 1.0) * 2π)

    δψ, _ = nutation(iau1980, ep; scale=scale)
    ϵ₀ = obliquity(iau1980, ep; scale=scale)

    return δψ * cos(ϵ₀) + sec2rad(0.00264 * sin(ω) + 0.000063 * sin(ω + ω))
end

function equinoxes(::IAU2000, ep::Epoch; scale=TT)
    t = julian_period(ep; scale=scale, unit=centuries, raw=true)

    l = fundamental(luna, t)
    lp = fundamental(sun, t)
    f = fundamental(luna, Longitude(), t)
    d = fundamental(luna, Elongation(), t)
    om = fundamental(luna, AscendingNode(), t)
    ve = fundamental(venus, t)
    le = fundamental(earth, t)
    pa = fundamental(t)

    s0 = 0.0
    s1 = 0.0

    for x in reverse(E0)
        eq = (
            x.nl * l +
            x.nlp * lp +
            x.nf * f +
            x.nd * d +
            x.nom * om +
            x.nve * ve +
            x.nle * le +
            x.npa * pa
        )
        se, ce = sincos(eq)
        s0 += x.s * se + x.c * ce
    end

    for x in reverse(E1)
        eq = (
            x.nl * l +
            x.nlp * lp +
            x.nf * f +
            x.nd * d +
            x.nom * om +
            x.nve * ve +
            x.nle * le +
            x.npa * pa
        )
        se, ce = sincos(eq)
        s1 += x.s * se + x.c * ce
    end

    return sec2rad(s0 + s1 * t)
end

function equinoxes(::IAU2000, ep::Epoch, ϵ, δψ; scale=TT)
    δψ * cos(ϵ) + equinoxes(iau2000, ep; scale=scale)
end

function equinoxes(iau::IAU2000Model, ep::Epoch; scale=TT)
    _, δϵ_pr = precession(iau2000, ep; scale=scale)
    ϵ = obliquity(iau1980, ep; scale=scale) + δϵ_pr
    δψ, _ = nutation(iau, ep; scale=scale)
    return equinoxes(iau2000, ep, ϵ, δψ; scale=scale)
end

function equinoxes(::IAU2006A, ep::Epoch)
    gst06a = apparent_sidereal(iau2006a, ep)
    gmst06 = mean_sidereal(iau2006, ep)
    return normalize_angle(gst06a - gmst06)
end

const SECS_TO_RAD = 7.272205216643039903848712e-5

function mean_sidereal(::IAU1982, ep::Epoch)
    t = julian_period(ep; scale=UT1, unit=centuries, raw=true)
    f = julian_period(ep; scale=UT1, raw=true) % 1.0 * SECONDS_PER_DAY
    gmst0 = @evalpoly(
        t,
        24110.54841 - SECONDS_PER_DAY / 2.0 + f,
        8640184.812866,
        0.093104,
        -6.2e-6,
    )
    return mod2pi(SECS_TO_RAD * gmst0)
end

function mean_sidereal(::IAU2000, ep::Epoch; scale=TT)
    t = julian_period(ep; scale=scale, unit=centuries, raw=true)
    gmst0 = @evalpoly(
        t,
        0.014506,
        4612.15739966,
        1.39667721,
        -0.00009344,
        0.00001882,
    )
    return mod2pi(earth_rotation_angle(iau2000, ep) + sec2rad(gmst0))
end

function mean_sidereal(::IAU2006, ep::Epoch)
    t = ep |> TTEpoch |> j2000 |> centuries |> value
    gmst0 = @evalpoly(
        t,
        0.014506,
        4612.156534,
        1.3915817,
        -0.00000044,
        -0.000029956,
        -0.0000000368,
    )
    return mod2pi(earth_rotation_angle(iau2000, ep) + sec2rad(gmst0))
end

function apparent_sidereal(::IAU1994, ep::Epoch)
    gmst82 = mean_sidereal(iau1982, ep)
    eqeq94 = equinoxes(iau1994, ep; scale=UT1)
    return mod2pi(gmst82 + eqeq94)
end

function apparent_sidereal(::IAU2000A, ep::Epoch)
    ee = equinoxes(iau2000a, ep)
    gmst = mean_sidereal(iau2000, ep)
    return mod2pi(gmst + ee)
end

function apparent_sidereal(::IAU2000B, ep::Epoch)
    ee = equinoxes(iau2000b, ep; scale=UT1)
    gmst = mean_sidereal(iau2000, ep; scale=UT1)
    return mod2pi(gmst + ee)
end

function apparent_sidereal(::IAU2006, ep::Epoch, rnpb)
    x, y = cip_coords(rnpb)
    s = cio_locator(iau2006, ep, x, y)
    era = earth_rotation_angle(iau2000, ep)
    eors = equation_of_origins(rnpb, s)
    return mod2pi(era - eors)
end

function apparent_sidereal(::IAU2006A, ep::Epoch)
    rnpb = precession_nutation_matrix(iau2006a, ep)
    return apparent_sidereal(iau2006, ep, rnpb)
end

