using AstroTime: Epoch, TTEpoch, centuries, j2000, value

export nutation

include(joinpath("constants", "nutation_2000b.jl"))

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

        dp += (x.ps + x.pst * t) * sarg + x.pc * carg
        de += (x.ec + x.ect * t) * carg + x.es * sarg
    end

    # Luni-solar terms
    δψ_ls = sec2rad(dp * 1e-7)
    δϵ_ls = sec2rad(de * 1e-7)

    # Fixed offsets in lieu of planetary terms
    δψ_pl = sec2rad(-0.135 * 1e-3)
    δϵ_pl = sec2rad(0.388 * 1e-3)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end
