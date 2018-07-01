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
const J000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const TURNAS = 1296000.0
const ARC_TO_RAD = 4.848136811095359935899141e-6

fal03  = t -> mod(485868.249036 + t * (1717915923.2178 + t * (31.8792 +t * (0.051635 +t * (- 0.00024470 )))), TURNAS) * ARC_TO_RAD
falp03 = t -> mod(1287104.793048 + t * ( 129596581.0481 + t * (- 0.5532 + t * (0.000136 + t * (- 0.00001149)))), TURNAS) * ARC_TO_RAD
faf03  = t -> mod(335779.526232 + t * ( 1739527262.8478 + t * (- 12.7512 + t * (- 0.001037 + t * (0.00000417)))), TURNAS) * ARC_TO_RAD
fad03  = t -> mod(1072260.703692 + t * ( 1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169 )))), TURNAS) * ARC_TO_RAD
faom03 = t -> mod(450160.398036 + t * ( -6962890.5431 + t * ( 7.4722 + t * ( 0.007702 + t * ( -0.00005939 )))), TURNAS ) * ARC_TO_RAD
fame03 = t -> mod(4.402608842 + 2608.7903141574 * t, 2 * pi)
fave03 = t -> mod(3.176146697 + 1021.3285546211 * t, 2 * pi)
fae03  = t -> mod(1.753470314 + 628.3075849991 * t, 2 * pi)
fama03 = t -> mod(6.203480913 + 334.0612426700 * t, 2 * pi)
faju03 = t -> mod(0.599546497 + 52.9690962641 * t, 2 * pi)
fasa03 = t -> mod(0.874016757 + 21.3299104960 * t, 2 * pi)
faur03 = t -> mod(5.481293872 + 7.4781598567 * t, 2 * pi)
fane03 = t -> mod(5.311886287 + 3.8133035638 * t, 2 * pi)
fapa03 = t -> (0.024381750 + 0.00000538691 * t) * t

# function xy06()
#     t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC
#
#     # Powers of T.
#     w = 1.0
#     for (jpt = 0 jpt <= MAXPT jpt++) {
#     pt[jpt] = w
#     w *= t
#     }
#
#     # Initialize totals in X and Y:  polynomial, luni-solar, planetary.
#     for (jxy = 0 jxy < 2 jxy++) {
#     xypr[jxy] = 0.0
#     xyls[jxy] = 0.0
#     xypl[jxy] = 0.0
#     }
#
#     # ---------------------------------
#     # Fundamental arguments (IERS 2003)
#     # ---------------------------------
#
#     # Mean anomaly of the Moon.
#     fa[0] = fal03(t)
#
#     # Mean anomaly of the Sun.
#     fa[1] = falp03(t)
#
#     # Mean argument of the latitude of the Moon.
#     fa[2] = faf03(t)
#
#     # Mean elongation of the Moon from the Sun.
#     fa[3] = fad03(t)
#
#     # Mean longitude of the ascending node of the Moon.
#     fa[4] = faom03(t)
#
#     # Planetary longitudes, Mercury through Neptune.
#     fa[5] = fame03(t)
#     fa[6] = fave03(t)
#     fa[7] = fae03(t)
#     fa[8] = fama03(t)
#     fa[9] = faju03(t)
#     fa[10] = fasa03(t)
#     fa[11] = faur03(t)
#     fa[12] = fane03(t)
#
#     # General accumulated precession in longitude.
#     fa[13] = fapa03(t)
#
#     # --------------------------------------
#     # Polynomial part of precession-nutation
#     # --------------------------------------
#
#     for (jxy = 0 jxy < 2 jxy++) {
#     for (j = MAXPT j >= 0 j--) {
#        xypr[jxy] += xyp[jxy][j] * pt[j]
#     }
#     }
#
#     # ----------------------------------
#     # Nutation periodic terms, planetary
#     # ----------------------------------
#
#     # Work backwards through the coefficients per frequency list.
#     ialast = NA
#     for (ifreq = NFPL-1 ifreq >= 0 ifreq--) {
#
#     # Obtain the argument functions.
#     arg = 0.0
#     for (i = 0 i < 14 i++) {
#        m = mfapl[ifreq][i]
#        if (m != 0) arg += (double)m * fa[i]
#     }
#     sc[0] = sin(arg)
#     sc[1] = cos(arg)
#
#     # Work backwards through the amplitudes at this frequency.
#     ia = nc[ifreq+NFLS]
#     for (i = ialast i >= ia i--) {
#
#     # Coefficient number (0 = 1st).
#        j = i-ia
#
#     # X or Y.
#        jxy = jaxy[j]
#
#     # Sin or cos.
#        jsc = jasc[j]
#
#     # Power of T.
#        jpt = japt[j]
#
#     # Accumulate the component.
#        xypl[jxy] += a[i-1] * sc[jsc] * pt[jpt]
#     }
#     ialast = ia-1
#     }
#
#     # -----------------------------------
#     # Nutation periodic terms, luni-solar
#     # -----------------------------------
#
#     # Continue working backwards through the number of coefficients list.
#     for (ifreq = NFLS-1 ifreq >= 0 ifreq--) {
#
#     # Obtain the argument functions.
#     arg = 0.0
#     for (i = 0 i < 5 i++) {
#        m = mfals[ifreq][i]
#        if (m != 0) arg += (double)m * fa[i]
#     }
#     sc[0] = sin(arg)
#     sc[1] = cos(arg)
#
#     # Work backwards through the amplitudes at this frequency.
#     ia = nc[ifreq]
#     for (i = ialast i >= ia i--) {
#
#     # Coefficient number (0 = 1st).
#        j = i-ia
#
#     # X or Y.
#        jxy = jaxy[j]
#
#     # Sin or cos.
#        jsc = jasc[j]
#
#     # Power of T.
#        jpt = japt[j]
#
#     # Accumulate the component.
#        xyls[jxy] += a[i-1] * sc[jsc] * pt[jpt]
#     }
#     ialast = ia-1
#     }
#
#     # ------------------------------------
#     # Results:  CIP unit vector components
#     # ------------------------------------
#
#     *x = ARC_TO_RAD * (xypr[0] + (xyls[0] + xypl[0]) / 1e6)
#     *y = ARC_TO_RAD * (xypr[1] + (xyls[1] + xypl[1]) / 1e6)
#
#
# end

end # module
