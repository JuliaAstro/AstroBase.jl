using ..Bodies

export fundamental

function fundamental(::Sun, t)
    p = @evalpoly(t, 1287104.793048, 129596581.0481, -0.5532, 0.000136, -0.00001149)
    return sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

function fundamental(::Luna, t)
    p = @evalpoly(t, 485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470)
    return sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

struct Longitude end

function fundamental(::Luna, ::Longitude, t)
    p = @evalpoly(t, 335779.526232, 1739527262.8478, -12.7512, -0.001037, 0.00000417)
    return sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

struct Elongation end

function fundamental(::Luna, ::Elongation, t)
    p = @evalpoly(t, 1072260.703692, 1602961601.2090, -6.3706, 0.006593, -0.00003169)
    return sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

struct AscendingNode end

function fundamental(::Luna, ::AscendingNode, t)
    p = @evalpoly(t, 450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939)
    return sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

fundamental(::Mercury, t) = mod2pi(4.402608842 + 2608.7903141574t)
fundamental(::Venus, t) = mod2pi(3.176146697 + 1021.3285546211t)
fundamental(::Earth, t) = mod2pi(1.753470314 + 628.3075849991t)
fundamental(::Mars, t) = mod2pi(6.203480913 + 334.0612426700t)
fundamental(::Jupiter, t) = mod2pi(0.599546497 + 52.9690962641t)
fundamental(::Saturn, t) = mod2pi(0.874016757 + 21.3299104960t)
fundamental(::Uranus, t) = mod2pi(5.481293872 + 7.4781598567t)
fundamental(::Neptune, t) = mod2pi(5.311886287 + 3.8133035638t)
fundamental(t) = @evalpoly(t, 0.0, 0.024381750, 0.00000538691)
