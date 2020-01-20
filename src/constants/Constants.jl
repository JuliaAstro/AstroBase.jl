module Constants

export astronomical_unit

const AU = 1.49597870700e8 # km

astronomical_unit(::Type{Float64}) = AU
astronomical_unit() = astronomical_unit(Float64)

end

