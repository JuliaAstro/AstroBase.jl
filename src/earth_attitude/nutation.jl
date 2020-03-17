using AstroTime: Epoch, TTEpoch, centuries, j2000, value
using ..Util: normalize_angle

export nutation

struct NutationLuniSolar
    # Coefficients of l,l',F,D,Om
    nl::Int
    nlp::Int
    nf::Int
    nd::Int
    nom::Int
    # longitude sin, t*sin, cos coefficients
    sp::Float64
    spt::Float64
    cp::Float64
    # obliquity cos, t*cos, sin coefficients
    ce::Float64
    cet::Float64
    se::Float64
end

struct NutationPlanetary
    # Coefficients of l, F, D and Omega
    nl::Int
    nf::Int
    nd::Int
    nom::Int
    # Coefficients of planetary longitudes
    nme::Int
    nve::Int
    nea::Int
    nma::Int
    nju::Int
    nsa::Int
    nur::Int
    nne::Int
    # Coefficient of general precession
    npa::Int
    # Longitude sin, cos coefficients
    sp::Int
    cp::Int
    # Obliquity sin, cos coefficients
    se::Int
    ce::Int
end

include(joinpath("constants", "nutation_1980.jl"))
include(joinpath("constants", "nutation_2000a.jl"))
include(joinpath("constants", "nutation_2000b.jl"))

"""
    nutation(jd1, jd2)

Returns nutation in longitude(radians) and obliquity(radians) for a given 2 part Julian date (TT format).

# Example

julia> nutation(2.4578265e6, 0.30434616919175345)
(-3.7565297299394694e-5, -3.665617105048724e-5
```
"""
function nutation(::IAU1980, ep)
    t = ep |> TTEpoch |> j2000 |> centuries |> value

    el_poly = @evalpoly t 485866.733 715922.633 31.310 0.064
    el = normalize_angle(sec2rad(el_poly) + (1325.0 * t % 1.0) * 2π)
    elp_poly = @evalpoly t 1287099.804 1292581.224 -0.577 -0.012
    elp = normalize_angle(sec2rad(elp_poly) + (99.0 * t % 1.0) * 2π)
    f_poly = @evalpoly t 335778.877 295263.137 -13.257 0.011
    f = normalize_angle(sec2rad(f_poly) + (1342.0 * t % 1.0) * 2π)
    d_poly = @evalpoly t 1072261.307 1105601.328 -6.891 0.019
    d = normalize_angle(sec2rad(d_poly) + (1236.0 * t % 1.0) * 2π)
    om_poly = @evalpoly t 450160.280 -482890.539 7.455 0.008
    om = normalize_angle(sec2rad(om_poly) + (-5.0 * t % 1.0) * 2π)

    dp = 0.0
    de = 0.0

    for x in reverse(NUTATION_1980)
        arg = x.nl * el + x.nlp * elp + x.nf  * f + x.nd  * d + x.nom * om
        sarg, carg = sincos(arg)

        s = x.sp + x.spt * t
        c = x.ce + x.cet * t
        dp += s * sarg
        de += c * carg
    end

    δψ = sec2rad(dp * 1e-4)
    δϵ = sec2rad(de * 1e-4)

    return δψ, δϵ
end

"""
    nutation_00a(jd1, jd2)

Returns luni-solar and planetary nutations for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> nutation_00a(2.4578265e6, 0.30440190993249416)
(-3.7589903391357206e-5, -3.6659049617818334e-5)
```
"""
function nutation(::IAU2000A, ep::Epoch)
    t = ep |> TTEpoch |> j2000 |> centuries |> value

    el = mean_anomaly(luna, t)
    f  = mean_longitude_minus_lan(luna, t)
    om = mean_longitude_ascending_node(luna, t)

    elp_as = @evalpoly t 1287104.79305 129596581.0481 -0.5532 0.000136 -0.00001149
    elp = sec2rad(elp_as % ARCSECONDS_IN_CIRCLE)
    d_as = @evalpoly t 1072260.70369 1602961601.2090 -6.3706 0.006593 -0.00003169
    d   = sec2rad(d_as % ARCSECONDS_IN_CIRCLE)

    dpls = 0.0
    dels = 0.0

    for x in reverse(NUTATION_2000A_LS)
        arg = mod2pi(x.nl * el + x.nlp * elp + x.nf * f + x.nd * d + x.nom * om)
        sarg, carg = sincos(arg)

        dpls += (x.sp + x.spt * t) * sarg + x.cp * carg
        dels += (x.ce + x.cet * t) * carg + x.se * sarg
    end

    al   = mod2pi(2.35555598 + 8328.6914269554 * t)
    af   = mod2pi(1.627905234 + 8433.466158131 * t)
    ad   = mod2pi(5.198466741 + 7771.3771468121 * t)
    aom  = mod2pi(2.18243920 - 33.757045 * t)
    alne = mod2pi(5.321159000 + 3.8127774000 * t)

    apa  = general_precession_in_longitude(t)
    alme = mean_longitude(mercury, t)
    alve = mean_longitude(venus, t)
    alea = mean_longitude(earth, t)
    alma = mean_longitude(mars, t)
    alju = mean_longitude(jupiter, t)
    alsa = mean_longitude(saturn, t)
    alur = mean_longitude(uranus, t)

    dppl = 0.0
    depl = 0.0

    for x in reverse(NUTATION_2000A_PL)
        arg = mod2pi(x.nl * al +
                     x.nf * af   +
                     x.nd * ad   +
                     x.nom * aom  +
                     x.nme * alme +
                     x.nve * alve +
                     x.nea * alea +
                     x.nma * alma +
                     x.nju * alju +
                     x.nsa * alsa +
                     x.nur * alur +
                     x.nne * alne +
                     x.npa * apa)
        sarg, carg = sincos(arg)

        dppl += x.sp * sarg + x.cp * carg
        depl += x.se * sarg + x.ce * carg
    end

    # Luni-solar terms
    δψ_ls = sec2rad(dpls * 1e-7)
    δϵ_ls = sec2rad(dels * 1e-7)

    # Planetary terms
    δψ_pl = sec2rad(dppl * 1e-7)
    δϵ_pl = sec2rad(depl * 1e-7)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end

"""
    nutation_00b(jd1, jd2)

Returns luni-solar and planetary nutation for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> nutation_00b(2.4578265e6, 0.30440190993249416)
(-3.7589177912131684e-5, -3.6657431214029895e-5)
```
"""
function nutation(::IAU2000B, ep::Epoch)
    t = ep |> TTEpoch |> j2000 |> centuries |> value

    el  = (485868.249036 + (1717915923.2178) * t) % ARCSECONDS_IN_CIRCLE |> sec2rad
    elp = (1287104.79305 + (129596581.0481) * t) % ARCSECONDS_IN_CIRCLE |> sec2rad
    f   = (335779.526232 + (1739527262.8478) * t) % ARCSECONDS_IN_CIRCLE |> sec2rad
    d   = (1072260.70369 + (1602961601.2090) * t) % ARCSECONDS_IN_CIRCLE |> sec2rad
    om  = (450160.398036 + (-6962890.5431) * t) % ARCSECONDS_IN_CIRCLE |> sec2rad

    dp = 0.0
    de = 0.0

    for x in reverse(NUTATION_2000B)
        arg = mod2pi(x.nl * el + x.nlp * elp + x.nf * f + x.nd * d + x.nom * om)
        sarg, carg = sincos(arg)

        dp += (x.sp + x.spt * t) * sarg + x.cp * carg
        de += (x.ce + x.cet * t) * carg + x.se * sarg
    end

    # Luni-solar terms
    δψ_ls = sec2rad(dp * 1e-7)
    δϵ_ls = sec2rad(de * 1e-7)

    # Fixed offsets in lieu of planetary terms
    δψ_pl = sec2rad(-0.135 * 1e-3)
    δϵ_pl = sec2rad(0.388 * 1e-3)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end

function nutation(::IAU2006, ep::Epoch)
    t = ep |> TTEpoch |> j2000 |> centuries |> value

    fj2 = -2.7774e-6 * t
    δψ, δϵ = nutation(iau2000a, ep)

    δψ += δψ * (0.4697e-6 + fj2)
    δϵ += δϵ * fj2

    return δψ, δϵ
end

