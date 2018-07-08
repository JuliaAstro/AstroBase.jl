module AstroBase

# package code goes here
export fal03, falp03, faf03, fad03, faom03, fame03, fave03, fae03, fama03, faju03, fasa03,
        faur03, fane03, fapa03

include("mfals.jl")
# NFLS = (int) (sizeof mfals / sizeof (int) / 5)
# # Maximum power of T in the polynomials for X and Y
# enum { MAXPT = 5 }
# # Number of frequencies:  planetary
# static const int NFPL = (int) (sizeof fal03
# mfapl / sizeof (int) / 14)
const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const TURNAS = 1296000.0
const ARC_TO_RAD = 4.848136811095359935899141e-6

fal03(t)  =  mod(485868.249036 + t * (1717915923.2178 + t * (31.8792 +t * (0.051635 +t * (-0.00024470 )))), TURNAS) * ARC_TO_RAD
falp03(t) =  mod(1287104.793048 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))), TURNAS) * ARC_TO_RAD
faf03(t)  =  mod(335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))), TURNAS) * ARC_TO_RAD
fad03(t)  =  mod(1072260.703692 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169 )))), TURNAS) * ARC_TO_RAD
faom03(t) =  mod(450160.398036 + t * (-6962890.5431 + t * ( 7.4722 + t * ( 0.007702 + t * ( -0.00005939 )))), TURNAS ) * ARC_TO_RAD
fame03(t) =  mod2pi(4.402608842 + 2608.7903141574 * t)
fave03(t) =  mod2pi(3.176146697 + 1021.3285546211 * t)
fae03(t)  =  mod2pi(1.753470314 + 628.3075849991 * t)
fama03(t) =  mod2pi(6.203480913 + 334.0612426700 * t)
faju03(t) =  mod2pi(0.599546497 + 52.9690962641 * t)
fasa03(t) =  mod2pi(0.874016757 + 21.3299104960 * t)
faur03(t) =  mod2pi(5.481293872 + 7.4781598567 * t)
fane03(t) =  mod2pi(5.311886287 + 3.8133035638 * t)
fapa03(t) =  (0.024381750 + 0.00000538691 * t) * t

function xy06(jd1, jd2)
    NFLS, NFPL, NA, MAXPT= length(mfals), length(mfapl), length(a), 5
    pt, fa = Vector{Float64}(MAXPT), Vector{Float64}(14)
    xypr, xypl, xyls, sc = Vector{Float64}(2), Vector{Float64}(2), Vector{Float64}(2),Vector{Float64}(2)
    # Powers of T.
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    w = 1.0
    for i in range(1, MAXPT)
        pt[i] = w
        w *= t
    end

    # Initialize totals in X and Y:  polynomial, luni-solar, planetary.
    for i in range(1, 2)
        xypr[i] = 0.0
        xyls[i] = 0.0
        xypl[i] = 0.0
    end

    # ---------------------------------
    # Fundamental arguments (IERS 2003)
    # ---------------------------------

    # Mean anomaly of the Moon.
    fa[1] = fal03(t)

    # Mean anomaly of the Sun.
    fa[2] = falp03(t)

    # Mean argument of the latitude of the Moon.
    fa[3] = faf03(t)

    # Mean elongation of the Moon from the Sun.
    fa[4] = fad03(t)

    # Mean longitude of the ascending node of the Moon.
    fa[5] = faom03(t)

    # Planetary longitudes, Mercury through Neptune.
    fa[6]  = fame03(t)
    fa[7]  = fave03(t)
    fa[8]  = fae03(t)
    fa[9]  = fama03(t)
    fa[10]  = faju03(t)
    fa[11] = fasa03(t)
    fa[12] = faur03(t)
    fa[13] = fane03(t)

    # General accumulated precession in longitude.
    fa[14] = fapa03(t)

    # --------------------------------------
    # Polynomial part of precession-nutation
    # --------------------------------------

    for i in range(1,2)
        for j in reverse(range(1, MAXPT))
           xypr[i] += xyp[i][j] * pt[j]
        end
    end

    # ----------------------------------
    # Nutation periodic terms, planetary
    # ----------------------------------

    # Work backwards through the coefficients per frequency list.
    ialast = NA
    for ifreq in reverse(range(1, NFPL))
    # Obtain the argument functions.
        arg = 0.0
        for i in range(1,14)
           m = mfapl[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end

        sc[1] = sin(arg)
        sc[2] = cos(arg)

        # Work backwards through the amplitudes at this frequency.
        ia = nc[ifreq+NFLS]
        for i in reverse(range(ia, ia - ialast))

            # Coefficient number (0 = 1st).
               j = i - ia - 1

            # X or Y.
               jxy = jaxy[j]

            # Sin or cos.
               jsc = jasc[j]

            # Power of T.
               jpt = japt[j]

            # Accumulate the component.
               xypl[jxy] += a[i] * sc[jsc] * pt[jpt]
        end
        ialast = ia
    end

    # -----------------------------------
    # Nutation periodic terms, luni-solar
    # -----------------------------------

    # Continue working backwards through the number of coefficients list.
    for ifreq in reverse(range(1, NFLS))

        # Obtain the argument functions.
        arg = 0.0
        for i in range(1, 5)
           m = mfals[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end
        sc[1] = sin(arg)
        sc[2] = cos(arg)

        # Work backwards through the amplitudes at this frequency.
        ia = nc[ifreq]
        for i in reverse(range(ia, ia - ialast))

            # Coefficient number (0 = 1st).
               j = i - ia - 1

            # X or Y.
               jxy = jaxy[j]

            # Sin or cos.
               jsc = jasc[j]

            # Power of T.
               jpt = japt[j]

            # Accumulate the component.
               xyls[jxy] += a[i] * sc[jsc] * pt[jpt]
        end
        ialast = ia
    end

    # ------------------------------------
    # Results:  CIP unit vector components
    # ------------------------------------

    ARC_TO_RAD * (xypr[1] + (xyls[1] + xypl[1]) / 1e6), ARC_TO_RAD * (xypr[2] + (xyls[2] + xypl[2]) / 1e6)
end
end # module
