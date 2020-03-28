using AstroTime: UT1Epoch, fractionofday, j2000, value

export earth_rotation_angle

function earth_rotation_angle(::IAU2000, ep)
    t = ep |> UT1Epoch |> j2000 |> value
    f = t % 1.0
    return mod2pi(2π * (f + 0.7790572732640 + 0.00273781191135448t))
end

function equation_of_equinoxes_complementary_terms(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    fa = (
        fundamental(luna, t),
        fundamental(sun, t),
        fundamental(luna, Longitude(), t),
        fundamental(luna, Elongation(), t),
        fundamental(luna, AscendingNode(), t),
        fundamental(venus, t),
        fundamental(earth, t),
        fundamental(t),
    )

    s0 = 0.0
    s1 = 0.0

    for i in reverse(eachindex(e0_coefficent))
        a = 0.0
        for j in 1:8
            a += e0_coefficent[i][j] * fa[j]
        end
        s0 += e0_arg[i][1] * sin(a) + e0_arg[i][2] * cos(a)
    end

    for i in reverse(eachindex(e1_coefficent))
        a = 0.0;
        for j in 1:8
            a += e1_coefficent[i][j] * fa[j]
        end
        s1 += e1_arg[i][1] * sin(a) + e1_arg[i][2] * cos(a)
    end

    sec2rad(s0 + s1 * t)
end

function equation_of_equinoxes_94(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    om = rem2pi(sec2rad(@evalpoly t 450160.280 -482890.539 7.455 0.008)
      + mod(-5.0 * t, 1.0) * 2pi, RoundNearest)

    dpsi, deps = nutation(iau1980, TTEpoch(jd1 * days, jd2 * days, origin=:julian))
    eps0 = obliquity(iau1980, TTEpoch(jd1 * days, jd2 * days, origin=:julian))

    dpsi*cos(eps0) + sec2rad((0.00264 * sin(om) + 0.000063 * sin(om + om)))
end

function equation_of_equinoxes_00(jd1, jd2, epsa, dpsi)
    dpsi * cos(epsa) + equation_of_equinoxes_complementary_terms(jd1, jd2)
end

function equation_of_equinoxes_a00(jd1, jd2)
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    dpsipr, depspr = precession(iau2000, ep)
    epsa = obliquity(iau1980, ep) + depspr
    dpsi, deps = nutation(iau2000a, ep)
    dpsi * cos(epsa) + equation_of_equinoxes_complementary_terms(jd1, jd2)
end

function equation_of_equinoxes_b00(jd1, jd2)
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    dpsipr, depspr = precession(iau2000, ep)
    epsa = obliquity(iau1980, ep) + depspr
    dpsi, deps = nutation(iau2000b, ep)
    dpsi * cos(epsa) + equation_of_equinoxes_complementary_terms(jd1, jd2)
end

function equation_of_equinoxes_a06(jd1, jd2)
    gst06a = greenwich_mean_sidereal_time_a06(0.0, 0.0, jd1, jd2)
    gmst06 = greenwich_mean_sidereal_time06(0.0, 0.0, jd1, jd2)

    rem2pi(gst06a - gmst06, RoundNearest)
end

function greenwich_mean_sidereal_time82(jd1, jd2)
    A = 24110.54841  -  SECONDS_PER_DAY / 2.0
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    if jd1 < jd2
        d1 = jd1
        d2 = jd2
    else
        d1 = jd2
        d2 = jd1
    end

    t = (d1 + (d2 - J2000)) / DAYS_PER_CENTURY
    f = SECONDS_PER_DAY * (mod(d1, 1.0) + mod(d2, 1.0))
    mod2pi(7.272205216643039903848712e-5 * (@evalpoly t A + f B C D))
end

function greenwich_mean_sidereal_time00(ut1, ut2, tt1, tt2)
    t = ((tt1 - J2000) + tt2) / DAYS_PER_CENTURY
    ut = UT1Epoch(ut1 * days, ut2 * days, origin=:julian)
    mod2pi(earth_rotation_angle(iau2000, ut) +
           sec2rad(@evalpoly t 0.014506 4612.15739966 1.39667721 -0.00009344 0.00001882))
end

function greenwich_mean_sidereal_time06(ut1, ut2, tt1, tt2)
    t = ((tt1 - J2000) + tt2) / DAYS_PER_CENTURY
    ut = UT1Epoch(ut1 * days, ut2 * days, origin=:julian)
    mod2pi(earth_rotation_angle(iau2000, ut) +
           sec2rad(@evalpoly t 0.014506 4612.156534 1.3915817 -0.00000044 -0.000029956 -0.0000000368))
end

function greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
    x, y = cip_coords(rnpb)
    s = s06(tta, ttb, x, y)
    ut = UT1Epoch(uta * days, utb * days, origin=:julian)
    era = earth_rotation_angle(iau2000, ut)
    eors = equation_of_origins(rnpb,s)
    mod2pi(era - eors)
end

function greenwich_mean_sidereal_time_a06(uta, utb, tta, ttb)
    ep = TTEpoch(tta * days, ttb * days, origin=:julian)
    rnpb = precession_nutation_matrix(iau2006a, ep)
    greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
end

function greenwich_apparent_sidereal_time94(uta, utb)
    gmst82 = greenwich_mean_sidereal_time82(uta, utb)
    eqeq94 = equation_of_equinoxes_94(uta, utb)
    mod2pi(gmst82 + eqeq94)
end

function greenwich_apparent_sidereal_time_a00(uta, utb, tta, ttb)
    gmst00 = greenwich_mean_sidereal_time00(uta, utb, tta, ttb)
    ee00a = equation_of_equinoxes_a00(tta, ttb)
    mod2pi(gmst00 + ee00a)
end

function greenwich_apparent_sidereal_time_b00(uta, utb)
    gmst00 = greenwich_mean_sidereal_time00(uta, utb, uta, utb)
    ee00b = equation_of_equinoxes_b00(uta, utb)
    mod2pi(gmst00 + ee00b)
end

function greenwich_apparent_sidereal_time_a06(uta, utb, tta, ttb)
    ep = TTEpoch(tta * days, ttb * days, origin=:julian)
    rnpb = precession_nutation_matrix(iau2006a, ep)
    greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
end

