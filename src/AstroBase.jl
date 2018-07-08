module AstroBase

# package code goes here
export xy06, fal03, falp03, faf03, fad03, faom03, fame03, fave03, fae03, fama03, faju03, fasa03,
       faur03, fane03, fapa03

include("mfals.jl")

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const TURNAS = 1296000.0

sec2rad(sec::Real) = deg2rad(sec/3600)

fal03(t)  =  mod((@evalpoly t 485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470) , TURNAS) * sec2rad(1)
falp03(t) =  mod((@evalpoly t 1287104.793048 129596581.0481 -0.5532 0.000136 -0.00001149),  TURNAS) * sec2rad(1)
faf03(t)  =  mod((@evalpoly t 335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417), TURNAS) * sec2rad(1)
fad03(t)  =  mod((@evalpoly t 1072260.703692 1602961601.2090 -6.3706 0.006593 -0.00003169), TURNAS) * sec2rad(1)
faom03(t) =  mod((@evalpoly t 450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939), TURNAS ) * sec2rad(1)
fame03(t) =  mod2pi(4.402608842 + 2608.7903141574t)
fave03(t) =  mod2pi(3.176146697 + 1021.3285546211t)
fae03(t)  =  mod2pi(1.753470314 + 628.3075849991t)
fama03(t) =  mod2pi(6.203480913 + 334.0612426700t)
faju03(t) =  mod2pi(0.599546497 + 52.9690962641t)
fasa03(t) =  mod2pi(0.874016757 + 21.3299104960t)
faur03(t) =  mod2pi(5.481293872 + 7.4781598567t)
fane03(t) =  mod2pi(5.311886287 + 3.8133035638t)
fapa03(t) =  (0.024381750 + 0.00000538691t) * t

function xy06(jd1, jd2)
    NFLS, NFPL, NA, MAXPT= length(mfals), length(mfapl), length(amp), 5
    pt, fa = Vector{Float64}(MAXPT + 1), Vector{Float64}(14)
    xypr, xypl, xyls, sc = zeros(2), zeros(2), zeros(2), zeros(2)

    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    # Powers of T.
    w = 1.0
    for i in 1:MAXPT+1
        pt[i] = w
        w *= t
    end
    # ---------------------------------
    # Fundamental arguments (IERS 2003)
    # ---------------------------------
    fa[1]  = fal03(t)
    fa[2]  = falp03(t)
    fa[3]  = faf03(t)
    fa[4]  = fad03(t)
    fa[5]  = faom03(t)
    fa[6]  = fame03(t)
    fa[7]  = fave03(t)
    fa[8]  = fae03(t)
    fa[9]  = fama03(t)
    fa[10] = faju03(t)
    fa[11] = fasa03(t)
    fa[12] = faur03(t)
    fa[13] = fane03(t)
    fa[14] = fapa03(t)

    # --------------------------------------
    # Polynomial part of precession-nutation
    # --------------------------------------

    for i in 1:2
           xypr[i] = @evalpoly t xyp[1][1] xyp[1][2] xyp[1][3] xyp[1][4] xyp[1][5] xyp[1][6]
    end

    # ----------------------------------
    # Nutation periodic terms, planetary
    # ----------------------------------
    ialast = NA
    for ifreq in NFPL:-1:1
        arg = 0.0
        for i in range(1,14)
           m = mfapl[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end

        sc[2], sc[1] = reim(cis(arg))

        ia = nc[ifreq + NFLS]
        for i in ialast + 1: -1: ia + 1
               j = i - ia
               jxy = jaxy[j]
               jsc = jasc[j]
               jpt = japt[j]
               xypl[jxy] += amp[i-1] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    # -----------------------------------
    # Nutation periodic terms, luni-solar
    # -----------------------------------
    for ifreq in NFLS: -1: 1
        arg = 0.0
        for i in 1:5
           m = mfals[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end

        sc[2], sc[1] = reim(cis(arg))

        ia = nc[ifreq]
        for i in ialast + 1: -1: ia + 1
               j = i - ia
               jxy = jaxy[j]
               jsc = jasc[j]
               jpt = japt[j]
               xyls[jxy] += amp[i-1] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    # ------------------------------------
    # Results:  CIP unit vector components
    # ------------------------------------

    sec2rad(1) * (xypr[1] + (xyls[1] + xypl[1]) / 1e6), sec2rad(1) * (xypr[2] + (xyls[2] + xypl[2]) / 1e6)
end
end # module
