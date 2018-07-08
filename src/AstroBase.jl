module AstroBase

# package code goes here
export xy06, fal03, falp03, faf03, fad03, faom03, fame03, fave03, fae03, fama03, faju03, fasa03,
       faur03, fane03, fapa03

include("mfals.jl")

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const TURNAS = 1296000.0

"""
    sec2rad(sec::Real)

Returns radians for given seconds.

# Example

```jldoctest
julia> sec2rad(3600 * 30)
0.5235987755982988
```
"""
sec2rad(sec::Real) = deg2rad(sec/3600)

"""
    fal03(t::Real)

Returns mean anomaly of Moon for a given Julian century.

# Example

```jldoctest
julia> fal03(23.0)
0.5891752616281019
```
"""
fal03(t::Real)  =  mod((@evalpoly t 485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470) , TURNAS) * sec2rad(1)

"""
    falp03(t::Real)

Returns mean anomaly of Sun for a given Julian century.

# Example

```jldoctest
julia> falp03(23.0)
5.857396217361825
```
"""
falp03(t::Real) =  mod((@evalpoly t 1287104.793048 129596581.0481 -0.5532 0.000136 -0.00001149),  TURNAS) * sec2rad(1)

"""
    faf03(t::Real)

Returns mean longitude of the Moon for a given Julian century.

# Example

```jldoctest
julia> faf03(23.0)
3.103138156410118
```
"""
faf03(t::Real)  =  mod((@evalpoly t 335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417), TURNAS) * sec2rad(1)

"""
    fad03(t::Real)

Returns mean elongation of the Moon from the Sun for given Julian century.

# Example

```jldoctest
julia> fad03(23.0)
2.8012040574296484
```
"""
fad03(t::Real)  =  mod((@evalpoly t 1072260.703692 1602961601.2090 -6.3706 0.006593 -0.00003169), TURNAS) * sec2rad(1)

"""
    faom03(t::Real)

Return fundamental argument for a given Julian century.

# Example
```jldoctest
julia> faom03(23.0)
4.904897783682109
```
"""
faom03(t::Real) =  mod((@evalpoly t 450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939), TURNAS ) * sec2rad(1)

"""
    fame03(t::Real)

Returns mean longitude of Mercury for a given Julian century.

# Example

```jldoctest
julia> fame03(23.0)
2.160150897150834
```
"""
fame03(t::Real) =  mod2pi(4.402608842 + 2608.7903141574t)

"""
    fave03(t::Real)

Returns mean longitude of Venus for a given Julian century.

# Example

```jldoctest
julia> fave03(23.0)
0.9030394378238363
```
"""
fave03(t::Real) =  mod2pi(3.176146697 + 1021.3285546211t)

"""
    fae03(t::Real)

Returns mean longitude of Earth for a given Julian century.

# Example

```jldoctest
julia> fae03(23.0)
1.501718780251826
```
"""
fae03(t::Real)  =  mod2pi(1.753470314 + 628.3075849991t)

"""
    fama03(t::Real)

Returns mean longitude of Mars for a given Julian century.

# Example

```jldoctest
julia> fama03(23.0)
5.276431642365657
```
"""
fama03(t::Real) =  mod2pi(6.203480913 + 334.0612426700t)

"""
    faju03(t::Real)

Returns mean longitude of Jupiter for a given Julian century.

# Example

```jldoctest
julia> faju03(23.0)
6.233996285639864
```
"""
faju03(t::Real) =  mod2pi(0.599546497 + 52.9690962641t)

"""
    fasa03(t::Real)

Returns mean longitude of Saturn for a given Julian century.

# Example

```jldoctest
julia> fasa03(23.0)
1.3735042049922535
```
"""
fasa03(t::Real) =  mod2pi(0.874016757 + 21.3299104960t)

"""
    faur03(t::Real)

Returns mean longitude of Uranus for a given Julian century.

# Example

```jldoctest
julia> faur03(23.0)
1.5497819750715893
```
"""
faur03(t::Real) =  mod2pi(5.481293872 + 7.4781598567t)

"""
    fane03(t::Real)

Returns mean longitude of Neptune for a given Julian century.

# Example

```jldoctest
julia> fane03(23.0)
5.053273953885775
```
"""
fane03(t::Real) =  mod2pi(5.311886287 + 3.8133035638t)

"""
    fapa03(t::Real)

Returns general accumulated precession in longitude for a given Julian century.

# Example

```jldoctest
julia> fapa03(23.0)
0.56362992539
```
"""
fapa03(t::Real) =  (0.024381750 + 0.00000538691t) * t

"""
    xy06(jd1, jd2)

Returns X, Y coordinates of celestial intermediate pole for a given 2-part Julian date.

# Example

```jldoctest
julia> xy06(2.4578265e6, 0.30440190993249416)
(0.0016558850230577835, -3.986943362456243e-5)
```
"""
function xy06(jd1, jd2)
    NFLS = length(mfals)
    NFPL = length(mfapl)
    NA = length(amp)
    MAXPT= 5
    pt = Vector{Float64}(MAXPT + 1)
    fa = Vector{Float64}(14)
    xypr, xypl, xyls, sc = zeros(2), zeros(2), zeros(2), zeros(2)

    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    # Powers of T.
    w = 1.0
    for i in 1:MAXPT+1
        pt[i] = w
        w *= t
    end

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

    for i in 1:2
           xypr[i] = @evalpoly t xyp[i][1] xyp[i][2] xyp[i][3] xyp[i][4] xyp[i][5] xyp[i][6]
    end

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
        for i in ialast+1:-1:ia+1
               j = i - ia
               jxy = jaxy[j]
               jsc = jasc[j]
               jpt = japt[j]
               xypl[jxy] += amp[i-1] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    for ifreq in NFLS:-1:1
        arg = 0.0
        for i in 1:5
           m = mfals[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end

        sc[2], sc[1] = reim(cis(arg))

        ia = nc[ifreq]
        for i in ialast+1:-1:ia+1
               j = i - ia
               jxy = jaxy[j]
               jsc = jasc[j]
               jpt = japt[j]
               xyls[jxy] += amp[i-1] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    sec2rad(1) * (xypr[1] + (xyls[1] + xypl[1]) / 1e6), sec2rad(1) * (xypr[2] + (xyls[2] + xypl[2]) / 1e6)
end
end # module
