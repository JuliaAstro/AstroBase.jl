using ReferenceFrameRotations: angle_to_dcm, angleaxis_to_dcm, compose_rotation

using ..Time: Epoch, TT, centuries, julian_period
using ..Util: sec2rad

export
    celestial_to_intermediate,
    celestial_to_terrestrial_cio,
    celestial_to_terrestrial_equinox,
    cio_locator,
    cip_coords,
    cip_coords_cio_locator,
    equation_of_origins,
    polar_motion,
    tio_locator

function equation_of_origins(rnpb, s)
    x = rnpb[3,1]
    ax = x / (1.0 + rnpb[3,3])
    xs = 1.0 - ax * x
    ys = -ax * rnpb[3,2]
    zs = -x
    p = rnpb[1,1] * xs + rnpb[1,2] * ys + rnpb[1,3] * zs
    q = rnpb[2,1] * xs + rnpb[2,2] * ys + rnpb[2,3] * zs
    return ifelse(p != 0 || q != 0, s - atan(q, p), s)
end

"""
    cip_coords(rbpn)

Extract from the bias-precession-nutation matrix the X,Y coordinates
of the Celestial Intermediate Pole.

# References

- [ERFA](https://github.com/liberfa/erfa/blob/master/src/bpn2xy.c)
"""
cip_coords(rbpn::AbstractMatrix) = rbpn[3,1], rbpn[3,2]

include(joinpath("constants", "cip_coords.jl"))

function cip_coords(::IAU2006, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)

    # Powers of T.
    pt = [t^i for i = 0:6]

    fa = (fundamental(luna, t),
          fundamental(sun, t),
          fundamental(luna, Longitude(), t),
          fundamental(luna, Elongation(), t),
          fundamental(luna, AscendingNode(), t),
          fundamental(mercury, t),
          fundamental(venus, t),
          fundamental(earth, t),
          fundamental(mars, t),
          fundamental(jupiter, t),
          fundamental(saturn, t),
          fundamental(uranus, t),
          fundamental(neptune, t),
          fundamental(t))

    xpr = @evalpoly(
        t,
        -0.016617,
        2004.191898,
        -0.4297829,
        -0.19861834,
        0.000007578,
        0.0000059285,
    )

    ypr = @evalpoly(
        t,
        -0.006951,
        -0.025896,
        -22.4072747,
        0.00190059,
        0.001112526,
        0.0000001358,
    )

    xypl = zeros(2)
    xyls = zeros(2)

    ialast = NAMP
    for ifreq in reverse(eachindex(FUN_PL))
        arg = 0.0
        for i in eachindex(fa)
            m = FUN_PL[ifreq][i]
            arg += m * fa[i]
        end

        sc = sincos(arg)

        ia = POINTERS[ifreq + NLS]
        for i in ialast:-1:ia
            j = i - ia + 1
            jxy = JAXY[j]
            jsc = JASC[j]
            jpt = JAPT[j]
            xypl[jxy] += AMPLITUDES[i] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    for ifreq in reverse(eachindex(FUN_LS))
        arg = 0.0
        for i in 1:5
           m = FUN_LS[ifreq][i]
           arg += m * fa[i]
        end

        sc = sincos(arg)

        ia = POINTERS[ifreq]
        for i in ialast:-1:ia
            j = i - ia + 1
            jxy = JAXY[j]
            jsc = JASC[j]
            jpt = JAPT[j]
            xyls[jxy] += AMPLITUDES[i] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end
    x = sec2rad((xpr + (xyls[1] + xypl[1]) / 1e6))
    y = sec2rad(ypr + (xyls[2] + xypl[2]) / 1e6)

    return x, y
end

struct SeriesTerm
    # Coefficients of l,l',F,D,Om,LVe,LE,pA
    nl::Int
    nlp::Int
    nf::Int
    nd::Int
    nom::Int
    nve::Int
    nle::Int
    npa::Int
    # Sine and cosine coefficients
    s::Float64
    c::Float64
end

include(joinpath("constants", "cio_locator.jl"))

function cio_locator(coeffs, terms, ep, x, y)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)

    l = fundamental(luna, t)
    lp = fundamental(sun, t)
    f = fundamental(luna, Longitude(), t)
    d = fundamental(luna, Elongation(), t)
    om = fundamental(luna, AscendingNode(), t)
    ve = fundamental(venus, t)
    le = fundamental(earth, t)
    pa = fundamental(t)

    w = copy(coeffs)

    for (i, t) in enumerate(terms)
        for x in reverse(t)
            a = (
                x.nl * l +
                x.nlp * lp +
                x.nf * f +
                x.nd * d +
                x.nom * om +
                x.nve * ve +
                x.nle * le +
                x.npa * pa
            )
            sa, ca = sincos(a)
            w[i] += x.s * sa + x.c * ca
        end
    end

    w0 = @evalpoly(t, w[1], w[2], w[3], w[4], w[5], w[6])
    return sec2rad(w0) - x * y / 2.0
end

cio_locator(::IAU2000, ep::Epoch, x, y) = cio_locator(SP_00, S00, ep, x, y)
cio_locator(::IAU2006, ep::Epoch, x, y) = cio_locator(SP_06, S06, ep, x, y)

function cio_locator(iau::IAU2000Model, ep::Epoch)
    rbpn = precession_nutation_matrix(iau, ep)
    x, y = cip_coords(rbpn)
    return cio_locator(iau2000, ep, x, y)
end

function cio_locator(::IAU2006A, ep::Epoch)
    rbpn = precession_nutation_matrix(iau2006a, ep)
    x, y = cip_coords(rbpn)
    return cio_locator(iau2006, ep, x, y)
end

function cip_coords_cio_locator(iau::IAU2000Model, ep)
    rbpn = precession_nutation_matrix(iau, ep)
    x, y = cip_coords(rbpn)
    return x, y, cio_locator(iau2000, ep, x, y)
end

function tio_locator(::IAU2000, ep)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)
    return sec2rad(-47e-6 * t)
end

function polar_motion(::IAU2000, xp, yp, sp)
    return angle_to_dcm(sp, -xp, -yp, :ZYX)
end

function celestial_to_intermediate(x, y, s)
    r2 = x^2 + y^2
    e = ifelse(r2 > 0.0, atan(y, x), 0.0)
    d = atan(sqrt(r2 / (1.0 - r2)))
    return angle_to_dcm(e, d, -(e+s), :ZYZ)
end

function celestial_to_intermediate(::IAU2000, ep::Epoch, x, y)
    s = cio_locator(iau2000, ep, x, y)
    return celestial_to_intermediate(x, y, s)
end

function celestial_to_intermediate(::IAU2000, ep::Epoch, rbpn)
    x, y = cip_coords(rbpn)
    return celestial_to_intermediate(iau2000, ep, x, y)
end

function celestial_to_intermediate(iau::IAU2000Model, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)
    rbpn = precession_nutation_matrix(iau, ep)
    return celestial_to_intermediate(iau2000, ep, rbpn)
end

function celestial_to_terrestrial_cio(rc2i, era, rpom)
    rera = angleaxis_to_dcm(era, [0.0, 0.0, 1.0])
    return compose_rotation(rc2i, rera, rpom)
end

function celestial_to_terrestrial_equinox(rbpn, gst, rpom)
    rgst = angleaxis_to_dcm(gst, [0.0, 0.0, 1.0])
    return compose_rotation(rbpn, rgst, rpom)
end

function celestial_to_terrestrial_cio(::IAU2000, ep::Epoch, x, y, xp, yp)
    rc2i = celestial_to_intermediate(iau2000, ep, x, y)
    era = earth_rotation_angle(iau2000, ep)
    sp = tio_locator(iau2000, ep)
    rpom = polar_motion(iau2000, xp, yp, sp)

    return celestial_to_terrestrial_cio(rc2i, era, rpom)
end

function celestial_to_terrestrial_cio(::IAU2000A, ep::Epoch, xp, yp)
    rc2i = celestial_to_intermediate(iau2000a, ep)
    era = earth_rotation_angle(iau2000, ep)
    sp = tio_locator(iau2000, ep)
    rpom = polar_motion(iau2000, xp, yp, sp)

    return celestial_to_terrestrial_cio(rc2i, era, rpom)
end

function celestial_to_terrestrial_cio(iau::IAU2000B, ep::Epoch, xp, yp)
    rc2i = celestial_to_intermediate(iau2000b, ep)
    era = earth_rotation_angle(iau2000, ep)
    rpom = polar_motion(iau2000, xp, yp, 0.0)

    return celestial_to_terrestrial_cio(rc2i, era, rpom)
end

function celestial_to_terrestrial_equinox(::IAU2000, ep::Epoch, δψ, δϵ, xp, yp)
    ϵ, _, _, _, _, rbpn = precession_nutation(iau2000, ep, δψ, δϵ)
    gmst = mean_sidereal(iau2000, ep)
    ee = equinoxes(iau2000, ep, ϵ, δψ)
    sp = tio_locator(iau2000, ep)
    rpom = polar_motion(iau2000, xp, yp, sp)

    return celestial_to_terrestrial_equinox(rbpn, gmst + ee, rpom)
end

